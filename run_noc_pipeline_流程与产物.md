# Phase 0 前后的通信数据流与输出格式

本文只说明 `run_noc_pipeline.sh` 中与通信描述有关的数据：

- DUMPI trace；
- CCDG；
- Phase 0 新增的 WSE plan；
- 它们之间的输入、输出和对齐关系。

## 1. Phase 0 改造前

```text
in.lammps
  → LAMMPS + DUMPI
  → dumpi-*.bin + dumpi-*.meta
  → dumpi2ccdg
  → raw / compact / trim-only CCDG
```

原流程只能从 MPI trace 恢复通信图。

## 2. Phase 0 改造后

```text
in.lammps
  → LAMMPS
       ├→ DUMPI trace → dumpi2ccdg → CCDG
       └→ 源码 emitter → rank JSONL → wse_plan.json

wse_plan.json + trim-only CCDG
  → validate_wse_plan.py
  → wse_plan_validation.json
```

改造后形成两条通信描述：

- CCDG：描述实际发生的 MPI 调用、计算间隔和依赖；
- WSE plan：描述 LAMMPS 算法直接产生的通信需求。

## 3. DUMPI 输出形式

每个 rank 一个二进制文件：

```text
dumpi-<timestamp>-0000.bin
dumpi-<timestamp>-0001.bin
...
dumpi-<timestamp>-0015.bin
```

内容是二进制 MPI 事件流，主要包括：

```text
MPI 操作类型
开始和结束时间
src / dst / tag / communicator
count / datatype
非阻塞 request ID
```

另有一个文本索引文件：

```text
dumpi-<timestamp>.meta
```

形式如下：

```text
hostname=<host>
numprocs=16
username=<user>
startime=<time>
fileprefix=dumpi-<timestamp>
version=16
subversion=0
subsubversion=0
```

`.meta` 给出 rank 数和 `.bin` 文件前缀；`dumpi2ccdg` 使用 `libundumpi`
读取这些文件。

## 4. CCDG 输出形式

`.ccdg` 文件是 JSON 文本，基本结构为：

```json
{
  "num_ranks": 16,
  "nodes": [],
  "cross_rank_edges": []
}
```

### 4.1 顶层字段

- `num_ranks`：MPI rank 总数；
- `nodes`：所有 rank 的计算和通信节点；
- `cross_rank_edges`：跨 rank 的发送—接收依赖。

### 4.2 计算节点

```json
{
  "id": 1,
  "rank": 0,
  "type": "COMPUTE",
  "compute_cycles": 9627,
  "compute_time_sec": 0.000003851,
  "compute_ops": 9627,
  "predecessors": [0]
}
```

含义：

- `id`：全局唯一节点编号；
- `rank`：节点属于哪个 rank；
- `compute_time_sec`：两次 MPI 调用之间的 CPU 时间；
- `compute_cycles`：按 CPU 频率换算的周期；
- `compute_ops`：当前模型使用的计算工作量；
- `predecessors`：执行该节点前必须完成的节点。

### 4.3 发送节点

```json
{
  "id": 6,
  "rank": 0,
  "type": "SEND",
  "comm_src": 0,
  "comm_dst": 12,
  "comm_tag": 0,
  "comm_bytes": 4128,
  "comm_count": 516,
  "predecessors": [5]
}
```

含义：

- rank 0 向 rank 12 发送；
- MPI tag 为 0；
- 发送 516 个元素；
- 总大小为 4128 bytes。

### 4.4 非阻塞接收节点

```json
{
  "id": 8,
  "rank": 0,
  "type": "IRECV",
  "comm_src": 4,
  "comm_dst": 0,
  "comm_bytes": 3936,
  "comm_count": 492,
  "pending_req_id": 96,
  "predecessors": [7]
}
```

`pending_req_id` 将 `IRECV` 与后续的 `WAIT/WAITALL` 关联。

### 4.5 等待节点

```json
{
  "id": 10,
  "rank": 0,
  "type": "WAIT",
  "pending_req_id": 96,
  "predecessors": [9, 8]
}
```

表示当前 rank 必须等待 request 96 对应的非阻塞操作完成。

### 4.6 集体通信节点

```json
{
  "id": 2,
  "rank": 0,
  "type": "ALLREDUCE",
  "comm_bytes": 64,
  "predecessors": [1]
}
```

`comm_bytes` 使用当前 CCDG 的 collective 总数据口径：

```text
num_ranks × count × datatype_bytes
```

### 4.7 跨 rank 依赖

```json
{
  "src_node": 6,
  "dst_node": 320
}
```

表示：

```text
发送节点 6
  → 远端接收或等待节点 320
```

BookSim 根据该边判断远端 rank 何时可以继续执行。

### 4.8 常见节点类型

```text
COMPUTE
SEND / ISEND
RECV / IRECV
WAIT / WAITALL
BARRIER
ALLREDUCE
BCAST
```

## 5. 三种 CCDG

### 5.1 Raw

```text
trace_<R>ranks_global.ccdg
```

形式：

```json
{
  "num_ranks": 16,
  "nodes": ["setup 节点", "run 节点"],
  "cross_rank_edges": ["全部依赖边"]
}
```

表示 setup 和 run 的完整 MPI 时间线。

### 5.2 Compact

```text
compact_<R>ranks_global.ccdg
```

形式：

```json
{
  "num_ranks": 16,
  "nodes": ["裁剪 setup 并折叠通信 burst 后的节点"],
  "cross_rank_edges": ["折叠后重新连接的依赖"]
}
```

连续通信可以合并为较大的节点：

```json
{
  "type": "SEND",
  "comm_src": 0,
  "comm_dst": 1,
  "comm_count": 9588,
  "comm_bytes": 76704
}
```

主要用于验证压缩前后通信字节守恒。

### 5.3 Trim-only

```text
trimonly_<R>ranks_global.ccdg
```

形式：

```json
{
  "num_ranks": 16,
  "nodes": ["只包含 run 段的原始粒度节点"],
  "cross_rank_edges": ["run 段原始依赖"]
}
```

它删除 setup，但不折叠消息，是 WSE plan 对拍和后续实验的主要载体。

三者关系：

```text
Raw       = setup + run 原始节点
Trim-only = run 原始节点
Compact   = run 节点 + 通信 burst 折叠
```

## 6. Phase 0 新增的 WSE rank 分片

### 6.1 CommBrick metadata

```json
{
  "kind": "metadata",
  "schema_version": 1,
  "rank": 0,
  "num_ranks": 16,
  "procgrid": [4, 4, 1],
  "myloc": [0, 0, 0]
}
```

表示 rank 数、处理器网格和当前 rank 的逻辑位置。

### 6.2 Swap setup

```json
{
  "kind": "swap_setup",
  "seq": 0,
  "epoch": 0,
  "timestep": 0,
  "rank": 0,
  "swap": 0,
  "dimension": 0,
  "direction": -1,
  "round": 0,
  "send_proc": 12,
  "recv_proc": 4,
  "pbc": true,
  "ghost_width": 10.3
}
```

表示 CommBrick 的静态通信轮次：

- `dimension`：X/Y/Z 维度；
- `direction`：负向或正向；
- `round`：第几轮 ghost 传播；
- `send_proc/recv_proc`：相邻 rank；
- `pbc`：是否跨周期边界；
- `ghost_width`：ghost 区域宽度。

### 6.3 CommBrick message

```json
{
  "kind": "message",
  "seq": 8,
  "scope": "run",
  "phase": "forward",
  "timestep": 1,
  "rank": 0,
  "src": 0,
  "dst": 12,
  "atom_count": 172,
  "value_count": 516,
  "datatype_bytes": 8,
  "bytes": 4128
}
```

`phase` 可以是：

```text
borders
forward
reverse
pair_forward
pair_reverse
```

消息字节数：

```text
bytes = value_count × datatype_bytes
```

### 6.4 Kspace/FFT message

```json
{
  "kind": "message",
  "component": "kspace",
  "seq": 4,
  "scope": "run",
  "phase": "fft_remap",
  "timestep": 1,
  "rank": 0,
  "src": 0,
  "dst": 5,
  "value_count": 108,
  "datatype_bytes": 8,
  "bytes": 864
}
```

`phase` 可以是：

```text
grid_reverse
grid_forward
fft_remap
```

分别表示：

- PPPM 电荷密度回收；
- PPPM 电场传播；
- FFT decomposition 转置。

### 6.5 Kspace collective

```json
{
  "kind": "collective",
  "component": "kspace",
  "seq": 100,
  "scope": "run",
  "phase": "kspace_reduce",
  "operation": "allreduce",
  "timestep": 1,
  "rank": 0,
  "value_count": 6,
  "datatype_bytes": 8,
  "bytes": 768
}
```

表示 PPPM energy/virial 全局归约。

## 7. 合并后的 `wse_plan.json`

输入：

```text
wse_plan.rankNNNN.jsonl
wse_plan.kspace.rankNNNN.jsonl
```

输出：

```json
{
  "schema_version": 1,
  "num_ranks": 16,
  "procgrid": [4, 4, 1],
  "rank_metadata": [],
  "component_metadata": [],
  "records": [],
  "summary": {
    "run": {
      "fft_remap": {
        "messages": 900,
        "collectives": 0,
        "bytes": 4805568
      }
    }
  }
}
```

字段含义：

- `rank_metadata`：CommBrick rank metadata；
- `component_metadata`：Kspace rank metadata；
- `records`：所有 setup/run 记录；
- `summary`：按 scope 和 phase 汇总消息、collective 和字节数。

数据流：

```text
每个 rank 的 CommBrick JSONL
每个 rank 的 Kspace JSONL
  → merge_wse_plan.py
  → 全局 wse_plan.json
```

## 8. WSE plan 与 CCDG 对拍形式

输入：

```text
wse_plan.json
trimonly_<R>ranks_global.ccdg
```

输出：

```text
wse_plan_validation.json
```

形式：

```json
{
  "passed": true,
  "num_ranks": 16,
  "plan_total_messages": 1220,
  "ccdg_total_messages": 1220,
  "plan_total_bytes": 7823808,
  "ccdg_total_bytes": 7823808,
  "directions": [
    {
      "direction": "-1,0",
      "plan_messages": 120,
      "ccdg_messages": 120,
      "plan_bytes": 1134144,
      "ccdg_bytes": 1134144,
      "relative_error": 0.0,
      "passed": true
    }
  ],
  "kspace_collectives_covered": true
}
```

对拍关系：

```text
WSE message
  → 按 src/dst 转换为二维方向
  → 与 CCDG SEND/ISEND 比较消息数和字节数

WSE collective
  → 按 operation 和 bytes
  → 检查 CCDG 中是否存在足够数量的对应 collective
```

## 9. 最终数据流

```text
                    ┌→ dumpi rank 二进制
LAMMPS 通信 ────────┤
                    └→ 源码 rank JSONL

dumpi rank 二进制
  → dumpi2ccdg
  → CCDG nodes + cross_rank_edges

源码 rank JSONL
  → merge_wse_plan.py
  → WSE records + summary

Trim-only CCDG + WSE plan
  → validate_wse_plan.py
  → 逐方向消息/字节对齐结果
```

核心区别：

```text
CCDG node
  = MPI 时间线中的计算、通信和依赖事件

WSE message
  = LAMMPS 算法直接给出的 src、dst、phase 和数据量

cross_rank_edge
  = SEND 与远端 RECV/WAIT 之间的执行依赖

swap_setup
  = CommBrick 的静态邻居和传播轮次
```
