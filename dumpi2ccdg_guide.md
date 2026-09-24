# MPI Trace 转 CCDG 转换指南

## 1. 概述

`dumpi2ccdg` 是一个用于将 DUMPI 格式的 MPI 通信 trace 转换为通信-计算依赖图（Communication-Computation Dependency Graph, CCDG）的工具。该工具使用 `libundumpi` 库解析二进制 trace 文件，通过回调机制将 MPI 事件转换为 CCDG 节点和边，最终输出 JSON 格式的依赖图文件。

### 1.1 工具位置

- **源码**: `/work1/jiangtao/lammps_trace/dumpi2ccdg/dumpi2ccdg.cpp`
- **可执行文件**: `/work1/jiangtao/lammps_trace/dumpi2ccdg/dumpi2ccdg`
- **依赖库**: `libundumpi`（DUMPI 解析库）

### 1.2 转换目标

| 输入 | 输出 |
|------|------|
| DUMPI 二进制 trace（`dumpi-*.bin`） | CCDG JSON 文件 |
| DUMPI 元数据（`dumpi-*.meta`） | 包含节点、边和依赖关系 |

---

## 2. 核心数据结构

### 2.1 CCDGNode（CCDG 节点）

每个节点代表一个计算或通信事件：

| 字段 | 类型 | 说明 |
|------|------|------|
| `id` | `uint64_t` | 全局唯一节点 ID |
| `rank` | `int` | 所属 MPI rank |
| `type` | `NodeType` | 节点类型（见下表） |
| `compute_cycles` | `double` | 计算周期（假设 2.5 GHz） |
| `compute_time_sec` | `double` | 计算时间（秒） |
| `comm_src` | `int` | 发送方 rank（通信节点） |
| `comm_dst` | `int` | 接收方 rank（通信节点） |
| `comm_tag` | `int` | MPI tag |
| `comm_bytes` | `uint64_t` | 通信字节数 |
| `collective_root` | `int` | 集体通信 root rank |
| `predecessors` | `vector<uint64_t>` | 同 rank 的前驱节点 ID |

### 2.2 NodeType（节点类型）

| 类型 | 对应 MPI 操作 |
|------|--------------|
| `COMPUTE` | MPI 调用之间的计算阶段 |
| `COMM_SEND` | `MPI_Send` |
| `COMM_RECV` | `MPI_Recv` |
| `COMM_ISEND` | `MPI_Isend` |
| `COMM_IRECV` | `MPI_Irecv` |
| `COMM_WAIT` | `MPI_Wait` |
| `COMM_WAITALL` | `MPI_Waitall` |
| `COMM_ALLREDUCE` | `MPI_Allreduce` |
| `COMM_BARRIER` | `MPI_Barrier` |
| `COMM_BCAST` | `MPI_Bcast` |
| `COMM_GATHER` | `MPI_Gather` |
| `COMM_ALLGATHER` | `MPI_Allgather` |
| `COMM_SCATTER` | `MPI_Scatter` |
| `COMM_ALLTOALL` | `MPI_Alltoall` |
| `COMM_REDUCE` | `MPI_Reduce` |

### 2.3 GlobalState（全局状态）

用于协调跨 rank 的依赖匹配：

| 字段 | 类型 | 说明 |
|------|------|------|
| `num_ranks` | `int` | MPI rank 数量 |
| `global_node_counter` | `uint64_t` | 全局节点计数器（确保 ID 唯一） |
| `ranks` | `vector<RankState>` | 每个 rank 的状态 |
| `pending_sends` | `map<tuple<int,int,int>, uint64_t>` | 待匹配的发送（src,dst,tag → node_id） |
| `ongoing_sends` | `map<pair<int,int>, uint64_t>` | 进行中的非阻塞发送 |
| `ongoing_recvs` | `map<pair<int,int>, uint64_t>` | 进行中的非阻塞接收 |
| `cross_edges` | `vector<pair<uint64_t,uint64_t>>` | 跨 rank 依赖边 |

---

## 3. 转换流程

转换分为四个阶段：

```
┌─────────────────────────────────────────────────────────────────────┐
│  阶段1: 查找 trace 文件                                              │
│  └─ find_trace_files() → 找到所有 dumpi-*.bin 和 dumpi-*.meta       │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│  阶段2: 解析每个 rank 的 trace                                        │
│  └─ parse_rank_trace() → 使用 libundumpi + 回调处理每个 MPI 调用     │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│  阶段3: 构建 CCDG 图                                                 │
│  ├─ COMPUTE 节点: MPI 调用之间的计算时间                              │
│  ├─ COMM 节点: Send/Recv/Isend/Irecv/Wait                            │
│  ├─ Collective 节点: Allreduce/Bcast/Barrier 等                      │
│  └─ 边: intra-rank (predecessors) + cross-rank (cross_rank_edges)   │
└──────────────────────────┬──────────────────────────────────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│  阶段4: 输出 JSON                                                    │
│  └─ print_json() → 输出 num_ranks, nodes[], cross_rank_edges[]      │
└─────────────────────────────────────────────────────────────────────┘
```

### 3.1 阶段1：查找 Trace 文件

函数 `find_trace_files()` 使用 `glob` 模式匹配：

- 查找 `dumpi-*.meta` 获取元数据文件
- 查找 `dumpi-*.bin` 获取各 rank 的二进制 trace 文件
- 按文件名排序确定 rank 顺序

### 3.2 阶段2：解析 Trace

函数 `parse_rank_trace()` 使用 `libundumpi` 库解析每个 rank 的 trace：

1. 打开 trace 文件：`undumpi_open(binfile)`
2. 读取 header：`undumpi_read_header(profile)`
3. 设置回调函数（见下表）
4. 读取 trace 流：`undumpi_read_stream(profile, &cb, &ctx, false)`
5. 读取 footer：`undumpi_read_footer(profile)`

### 3.3 阶段3：构建 CCDG

每个 MPI 调用触发对应的回调函数，构建节点和边：

#### MPI 事件回调映射

| MPI 操作 | 回调函数 | CCDG 节点类型 |
|----------|----------|--------------|
| `MPI_Send` | `on_send_cb()` | `COMM_SEND` |
| `MPI_Recv` | `on_recv_cb()` | `COMM_RECV` |
| `MPI_Isend` | `on_isend_cb()` | `COMM_ISEND` |
| `MPI_Irecv` | `on_irecv_cb()` | `COMM_IRECV` |
| `MPI_Wait` | `on_wait_cb()` | `COMM_WAIT` |
| `MPI_Waitall` | `on_waitall_cb()` | `COMM_WAITALL` |
| `MPI_Allreduce` | `on_allreduce_cb()` | `COMM_ALLREDUCE` |
| `MPI_Bcast` | `on_bcast_cb()` | `COMM_BCAST` |
| `MPI_Barrier` | `on_barrier_cb()` | `COMM_BARRIER` |
| `MPI_Gather` | `on_gather_cb()` | `COMM_GATHER` |
| `MPI_Allgather` | `on_allgather_cb()` | `COMM_ALLGATHER` |
| `MPI_Scatter` | `on_scatter_cb()` | `COMM_SCATTER` |
| `MPI_Alltoall` | `on_alltoall_cb()` | `COMM_ALLTOALL` |
| `MPI_Reduce` | `on_reduce_cb()` | `COMM_REDUCE` |

#### COMPUTE 节点生成逻辑

在每个 MPI 调用之前，如果存在时间间隔，则生成 COMPUTE 节点：

```
MPI_Send (t=10ms)  ─── 间隔 5ms ─── MPI_Recv (t=15ms)
     │                                      │
     │  add_compute_node(10ms, 15ms)        │
     │      ↓                               │
     └── COMPUTE(5ms) ──────────────────────┘
```

**时间计算**：
- `compute_time_sec = 当前 CPU 时间 - 上一次 MPI 调用的 CPU 时间`
- `compute_cycles = compute_time_sec × 2.5e9`（假设 CPU 频率为 2.5 GHz）

#### 点对点通信依赖匹配

**MPI_Send → MPI_Recv**：

```
Rank 0: MPI_Send(dest=1, tag=0)  ──→ 创建 SEND 节点
                                     存入 pending_sends[(0,1,0)] = node_id

Rank 1: MPI_Recv(source=0, tag=0) ──→ 创建 RECV 节点
                                     从 pending_sends[(0,1,0)] 查找
                                     添加跨 rank 边: SEND_node_id → RECV_node_id
```

**MPI_Isend → MPI_Wait**：

```
Rank 0: MPI_Isend(dest=1, tag=0) ──→ 创建 ISEND 节点
                                     记录到 pending_requests[req_id]

Rank 0: MPI_Wait(req)            ──→ 创建 WAIT 节点
                                     添加前驱: ISEND_node_id → WAIT_node_id
                                     如果是 Irecv 的 Wait，还需匹配远程 Send
```

#### 集体通信处理

集体通信直接创建对应类型的节点，不涉及跨 rank 边（由 BookSim 内部处理）：

```cpp
// MPI_Allreduce 示例
uint64_t bytes = prm->count * compute_datatype_size(prm->datatype);
// 集体通信的总字节数 = 单 rank 数据 × 总 rank 数
rs.add_collective_node(COMM_ALLREDUCE, t, bytes * num_ranks);
```

### 3.4 阶段4：输出 JSON

函数 `print_json()` 输出标准 JSON 格式：

```json
{
  "num_ranks": 4,
  "nodes": [...],
  "cross_rank_edges": [...]
}
```

---

## 4. 依赖边类型

### 4.1 Intra-rank 边（同 rank 内部依赖）

存储在每个节点的 `predecessors` 数组中，表示同 rank 内的执行顺序：

```
COMPUTE ──→ SEND ──→ COMPUTE ──→ RECV ──→ COMPUTE
  0           1         2          3         4
  ↓           ↓         ↓          ↓
[无]      [0]      [1]      [2]      [3]
```

### 4.2 Cross-rank 边（跨 rank 依赖）

存储在 `cross_rank_edges` 数组中，表示不同 rank 之间的通信依赖：

```
Rank 0: SEND node(5459)  ──cross_edge──→  Rank 1: RECV node(5737)
Rank 0: ISEND node(100)  ──cross_edge──→  Rank 1: WAIT node(200)
```

---

## 5. 编译与使用

### 5.1 编译

```bash
cd /work1/jiangtao/lammps_trace/dumpi2ccdg
make
```

### 5.2 使用方法

```bash
cd /work1/jiangtao/lammps_trace/dumpi2ccdg
./dumpi2ccdg <trace_directory> > <output.ccdg>
```

#### 示例1：转换 4-rank trace

```bash
./dumpi2ccdg ../runs/trace_4ranks_20260721_201258 > ../runs/trace_4ranks_20260721_201258/trace_4ranks_global.ccdg
```

#### 示例2：转换 16-rank trace

```bash
./dumpi2ccdg ../runs/trace_16ranks_20260721_201406 > ../runs/trace_16ranks_20260721_201406/trace_16ranks_global.ccdg
```

### 5.3 输入输出示例

**输入**：

```
runs/trace_4ranks_20260721_201258/
├── dumpi-2026.07.21.20.12.59-0000.bin  # Rank 0
├── dumpi-2026.07.21.20.12.59-0001.bin  # Rank 1
├── dumpi-2026.07.21.20.12.59-0002.bin  # Rank 2
├── dumpi-2026.07.21.20.12.59-0003.bin  # Rank 3
└── dumpi-2026.07.21.20.12.59.meta      # 元数据
```

**输出统计**：

```
Found 4 rank trace files in ../runs/trace_4ranks_20260721_201258/
Parsing rank 0/3: ../runs/trace_4ranks_20260721_201258/dumpi-2026.07.21.20.12.59-0000.bin
  -> 5570 nodes
Parsing rank 1/3: ../runs/trace_4ranks_20260721_201258/dumpi-2026.07.21.20.12.59-0001.bin
  -> 5538 nodes
Parsing rank 2/3: ../runs/trace_4ranks_20260721_201258/dumpi-2026.07.21.20.12.59-0002.bin
  -> 5538 nodes
Parsing rank 3/3: ../runs/trace_4ranks_20260721_201258/dumpi-2026.07.21.20.12.59-0003.bin
  -> 5538 nodes

Statistics:
  Total nodes: 22184
  Compute nodes: 11092
  Communication nodes: 11092
  Total compute time: 345.321 seconds
  Cross-rank edges: 841
```

---

## 6. CCDG JSON 输出格式

### 6.1 完整结构

```json
{
  "num_ranks": 4,
  "nodes": [
    {
      "id": 0,
      "rank": 0,
      "type": "COMPUTE",
      "compute_cycles": 90960475,
      "compute_time_sec": 0.03638419,
      "predecessors": []
    },
    {
      "id": 1,
      "rank": 0,
      "type": "BCAST",
      "comm_bytes": 16,
      "collective_root": 0,
      "predecessors": [0]
    },
    {
      "id": 5459,
      "rank": 0,
      "type": "SEND",
      "comm_src": 0,
      "comm_dst": 1,
      "comm_tag": 0,
      "comm_bytes": 16,
      "predecessors": [5458]
    },
    {
      "id": 5737,
      "rank": 1,
      "type": "RECV",
      "comm_src": 0,
      "comm_dst": 1,
      "comm_tag": 0,
      "comm_bytes": 16,
      "predecessors": [5736]
    }
  ],
  "cross_rank_edges": [
    {"src_node": 5459, "dst_node": 5737}
  ]
}
```

### 6.2 节点类型字段说明

| 节点类型 | 必选字段 | 可选字段 |
|----------|----------|----------|
| `COMPUTE` | `id`, `rank`, `type`, `compute_cycles` | `compute_time_sec`, `predecessors` |
| `SEND/RECV` | `id`, `rank`, `type`, `comm_src`, `comm_dst` | `comm_tag`, `comm_bytes`, `predecessors` |
| `ISEND/IRECV` | `id`, `rank`, `type`, `comm_src`, `comm_dst` | `comm_tag`, `comm_bytes`, `pending_req_id`, `predecessors` |
| `WAIT/WAITALL` | `id`, `rank`, `type` | `pending_req_id`, `predecessors` |
| 集体通信 | `id`, `rank`, `type`, `comm_bytes` | `collective_root`, `predecessors` |

---

## 7. 关键设计要点

| 设计点 | 实现方式 | 作用 |
|--------|----------|------|
| **全局节点 ID** | `global_node_counter` 统一递增 | 确保跨 rank 边引用唯一节点 |
| **依赖匹配** | `pending_sends[(src,dst,tag)]` | 匹配 Send-Recv 对 |
| **非阻塞请求跟踪** | `pending_requests[req_id]` | 关联 Isend/Irecv 和 Wait |
| **计算周期估算** | `cpu_time_sec × 2.5e9` | 将 CPU 时间转换为周期数 |
| **通信量计算** | `count × datatype_size` | 计算实际通信字节数 |
| **MPI 类型映射** | `compute_datatype_size(dumpi_datatype)` | 映射 MPI 数据类型到字节数 |

---

## 8. 与 BookSim 的集成

转换生成的 CCDG 文件直接用于 BookSim 2.0 的 `CCDGTrafficManager`：

1. **配置文件**：指定 CCDG 文件路径
   ```
   sim_type = ccdg;
   ccdg_file = /path/to/trace_4ranks_global.ccdg;
   ```

2. **仿真流程**：
   - `CCDGTrafficManager` 解析 CCDG JSON 文件
   - 为每个 rank 创建 PE（Processing Element）
   - PE 根据节点类型执行计算或注入网络数据包
   - 跨 rank 边用于同步不同 rank 的执行

---

## 9. 常见问题

### Q1: 为什么节点 ID 需要全局唯一？

**A**: 跨 rank 边需要引用具体的节点，如果每个 rank 独立编号会导致节点 ID 冲突，无法正确建立依赖关系。

### Q2: COMPUTE 节点的计算周期是如何估算的？

**A**: 根据 MPI 调用之间的 CPU 时间差计算，假设 CPU 频率为 2.5 GHz。实际频率应根据硬件环境调整。

### Q3: 非阻塞通信（Isend/Irecv）如何处理？

**A**: 使用 `pending_requests` 跟踪请求，当遇到 `Wait` 时，将对应的 ISEND/IRECV 节点设为前驱，并匹配远程的 Send 节点建立跨 rank 边。

### Q4: 集体通信为什么没有跨 rank 边？

**A**: 集体通信的同步语义在 BookSim 的 `CCDGTrafficManager` 内部处理，不需要显式的跨 rank 边。

### Q5: 如何查看原始 MPI 事件？

**A**: 使用 DUMPI 提供的 `dumpi2ascii` 工具：
```bash
cd /path/to/trace/directory
../../install/bin/dumpi2ascii -a -f dumpi-*-0000.bin | head -100
```

---

## 10. 参考文件

| 文件 | 路径 | 说明 |
|------|------|------|
| 转换工具源码 | `/work1/jiangtao/lammps_trace/dumpi2ccdg/dumpi2ccdg.cpp` | 核心转换逻辑 |
| BookSim 集成 | `/work1/jiangtao/lammps_trace/booksim2/src/ccdg_trafficmanager.cpp` | CCDG 仿真驱动 |
| 示例配置 | `/work1/jiangtao/lammps_trace/booksim2/ccdg_lammps_4x4.cfg` | BookSim 配置文件 |
| Trace 捕获脚本 | `/work1/jiangtao/lammps_trace/run_trace_capture.sh` | DUMPI trace 捕获 |