# Phase 0 前后 `run_noc_pipeline.sh` 各阶段产物

## 1. 总体变化

### Phase 0 改造前

```text
实验参数
  → Stage 1: in.lammps
  → Stage 2: LAMMPS + DUMPI → dumpi-*.bin/meta
  → Stage 3: dumpi2ccdg → raw/compact/trimonly CCDG
  → Stage 4: CCDG 质量检查
  → Stage 5: BookSim 五档实验
  → Stage 6: evaluation
```

通信需求只能从 DUMPI trace 间接恢复。

### Phase 0 改造后

```text
实验参数
  → Stage 1: in.lammps
  → Stage 2: LAMMPS
       ├→ DUMPI → dumpi-*.bin/meta
       └→ 源码 emitter → per-rank JSONL → wse_plan.json
  → Stage 3: dumpi2ccdg → raw/compact/trimonly CCDG
  → Stage 4: CCDG 质量检查 + wse_plan/CCDG 对拍
  → Stage 5: BookSim 五档实验
  → Stage 6: evaluation
```

Phase 0 主要改变 Stage 2 和 Stage 4。Stage 3、5、6 的原有产物格式保持不变，
因此现有实验仍可运行。

运行命令：

```bash
./run_noc_pipeline.sh <short|long> <ranks> <cu|h2o|lialocl> <atoms> [free|demand|both]
```

每次实验的产物目录：

```text
runs/pipeline/<mode>_<system>_<atoms>a_<ranks>r_<timestamp>/
```

下文简称 `RUN_DIR`。

## 2. Stage 1：生成 LAMMPS 输入

### 改造前产物

```text
RUN_DIR/
├── in.lammps
├── data.LiAlOCl_nvt_charge   # 仅 LiAlOCl
└── Cu_u3.eam                 # 仅 Cu
```

`in.lammps` 是文本格式，例如：

```text
units           metal
boundary        p p p
atom_style      charge
processors      4 4 1
pair_style      lj/cut/coul/long 10.0 10.0
kspace_style    pppm 1.0e-4
fix             1 all nve
run             1
```

含义：

- 定义原子体系、势函数和运行步数；
- `processors K K 1` 固定 Z 压缩处理器网格；
- short 使用截断库仑势；
- long 的 H2O/LiAlOCl 使用 PPPM Kspace/FFT；
- Cu 使用 EAM，不经过 Kspace。

### Phase 0 改造后产物

```text
与改造前相同
```

Phase 0 不改变实验输入，保证新旧通信描述来自完全相同的 LAMMPS 算例。

## 3. Stage 2：运行 LAMMPS 并捕获通信

## 3.1 改造前产物

### LAMMPS 日志

```text
RUN_DIR/
├── lammps.log
└── log.lammps
```

`lammps.log` 是终端输出文本：

```text
LAMMPS (...)
...
Loop time of 0.00184314 on 16 procs for 1 steps ...
```

含义：

- 判断 LAMMPS 是否成功；
- 提取正式 run 段的 `Loop time`；
- 发生错误时提供错误堆栈和退出信息。

`log.lammps` 是 LAMMPS 自身日志：

```text
4 by 4 by 1 MPI processor grid
...
Step Temp E_pair E_mol TotEng Press
```

含义：

- 检查实际处理器网格；
- 保存 thermo 和运行配置。

### DUMPI rank trace

```text
RUN_DIR/
├── dumpi-<timestamp>-0000.bin
├── dumpi-<timestamp>-0001.bin
├── ...
├── dumpi-<timestamp>-0015.bin
└── dumpi-<timestamp>.meta
```

`.bin` 是 DUMPI 二进制事件流，每个 rank 一个文件。内部记录：

```text
MPI 操作类型
开始/结束时间
src/dst/tag/communicator
count/datatype
非阻塞 request ID
```

它不能直接按普通文本阅读，需要 `libundumpi` 解析。

`.meta` 是文本索引：

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

含义：

- 指明 rank 总数和 `.bin` 文件前缀；
- 记录捕获环境和 DUMPI 格式版本；
- Stage 3 根据它找到所有 rank trace。

## 3.2 Phase 0 改造后新增产物

改造后保留以上全部文件，并新增源码级 WSE plan。

### CommBrick rank 分片

```text
RUN_DIR/
├── wse_plan.rank0000.jsonl
├── wse_plan.rank0001.jsonl
└── ...
```

格式为 JSON Lines，每行一个独立对象。

Metadata：

```json
{"kind":"metadata","schema_version":1,"rank":0,"num_ranks":16,
 "procgrid":[4,4,1],"myloc":[0,0,0]}
```

表示 rank 数、处理器网格和当前 rank 的逻辑坐标。

Swap setup：

```json
{"kind":"swap_setup","seq":0,"epoch":0,"timestep":0,"rank":0,
 "swap":0,"dimension":0,"direction":-1,"round":0,
 "send_proc":12,"recv_proc":4,"pbc":true,"ghost_width":10.3}
```

表示：

- 通信发生在哪个维度和方向；
- 属于第几轮 ghost 传播；
- 发送/接收邻居；
- 是否跨越 PBC；
- ghost 区域宽度。

Message：

```json
{"kind":"message","seq":8,"scope":"run","phase":"forward",
 "timestep":1,"rank":0,"src":0,"dst":12,
 "atom_count":172,"value_count":516,
 "datatype_bytes":8,"bytes":4128}
```

表示一条源码级点对点通信需求：

```text
bytes = value_count × datatype_bytes
```

`phase` 可以是：

```text
borders
forward
reverse
pair_forward
pair_reverse
```

### Kspace/FFT rank 分片

long 的 H2O/LiAlOCl 还会生成：

```text
RUN_DIR/
├── wse_plan.kspace.rank0000.jsonl
├── wse_plan.kspace.rank0001.jsonl
└── ...
```

点对点消息形式：

```json
{"kind":"message","component":"kspace","seq":4,
 "scope":"run","phase":"fft_remap","timestep":1,
 "rank":0,"src":0,"dst":5,
 "value_count":108,"datatype_bytes":8,"bytes":864}
```

`phase` 可以是：

```text
grid_reverse   # 电荷密度 ghost grid 向 owner 汇总
grid_forward   # 电场向 ghost grid 传播
fft_remap      # FFT decomposition 转置
```

Collective 形式：

```json
{"kind":"collective","component":"kspace","seq":100,
 "scope":"run","phase":"kspace_reduce","operation":"allreduce",
 "timestep":1,"rank":0,"value_count":6,
 "datatype_bytes":8,"bytes":768}
```

当前 collective 字节口径与 CCDG 一致：

```text
bytes = num_ranks × value_count × datatype_bytes
```

### 合并后的 `wse_plan.json`

```text
RUN_DIR/
├── wse_plan.json
└── wse_plan_merge.log
```

`wse_plan.json` 是普通 JSON：

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

含义：

- `rank_metadata`：CommBrick rank 信息；
- `component_metadata`：Kspace rank 信息；
- `records`：所有 setup/run 消息和 collective；
- `summary`：按 scope 和 phase 汇总消息数、collective 数和字节数。

`wse_plan_merge.log` 是文本：

```text
WSE plan merged: ranks=16 run_messages=1220
run_collectives=32 run_bytes=7838144 output=<path>/wse_plan.json
```

用于快速检查是否收齐全部 rank 分片。

## 4. Stage 3：DUMPI 转换成三种 CCDG

Phase 0 前后都生成相同格式的三种 CCDG。`.ccdg` 后缀对应的是 JSON 文本，
不是二进制文件。

## 4.1 CCDG 基本输出形式

```json
{
  "num_ranks": 16,
  "nodes": [
    {
      "id": 0,
      "rank": 0,
      "type": "BARRIER",
      "wall_time_sec": 15745836.766468445,
      "wall_duration_sec": 0.000002205
    },
    {
      "id": 1,
      "rank": 0,
      "type": "COMPUTE",
      "compute_cycles": 9627,
      "compute_time_sec": 0.000003851,
      "compute_ops": 9627,
      "predecessors": [0]
    },
    {
      "id": 2,
      "rank": 0,
      "type": "SEND",
      "comm_src": 0,
      "comm_dst": 12,
      "comm_tag": 0,
      "comm_bytes": 4128,
      "comm_count": 516,
      "predecessors": [1]
    }
  ],
  "cross_rank_edges": [
    {"src_node": 2, "dst_node": 320}
  ]
}
```

字段含义：

- `num_ranks`：参与图的 MPI rank 数；
- `nodes`：所有 rank 的计算和通信事件；
- `id`：全局唯一节点 ID；
- `rank`：节点所属 rank；
- `type`：COMPUTE、SEND、IRECV、WAIT、BARRIER、ALLREDUCE 等；
- `predecessors`：同 rank 内的执行依赖；
- `cross_rank_edges`：发送节点到远端接收/等待节点的数据依赖。

## 4.2 Raw CCDG

```text
RUN_DIR/
├── trace_<R>ranks_global.ccdg
└── ccdg_gen.log
```

输出形式：

```json
{
  "num_ranks": 16,
  "nodes": ["setup 节点", "run 节点"],
  "cross_rank_edges": ["全部跨 rank 边"]
}
```

含义：

- 保留 setup 和 run；
- 最接近完整 DUMPI 时间线；
- 文件最大；
- 用于回溯原始 MPI 行为，不直接作为主 BookSim 输入。

## 4.3 Compact CCDG

```text
RUN_DIR/
├── compact_<R>ranks_global.ccdg
└── ccdg_compact_gen.log
```

生成参数：

```text
CCDG_TRIM_SETUP=1
CCDG_COMPACT=1
```

输出仍是相同 JSON 结构：

```json
{
  "num_ranks": 16,
  "nodes": ["裁剪后并经过 burst 折叠的节点"],
  "cross_rank_edges": ["折叠后重新连接的依赖边"]
}
```

通信节点可能表示多条原始消息的合并结果：

```json
{
  "type": "SEND",
  "comm_src": 0,
  "comm_dst": 1,
  "comm_bytes": 76704,
  "comm_count": 9588
}
```

含义：

- 删除 setup；
- 合并连续通信 burst；
- `comm_bytes` 是合并后的总字节数；
- 用于检查 compaction 是否保持通信字节守恒；
- 当前不是五档 BookSim 的共同载体。

## 4.4 Trim-only CCDG

```text
RUN_DIR/
├── trimonly_<R>ranks_global.ccdg
└── ccdg_trimonly_gen.log
```

生成参数：

```text
CCDG_TRIM_SETUP=1
CCDG_COMPACT 未设置
```

输出形式：

```json
{
  "num_ranks": 16,
  "nodes": ["仅 run 段的原始粒度节点"],
  "cross_rank_edges": ["run 段的原始依赖边"]
}
```

含义：

- 删除 setup；
- 不折叠正式迭代中的 SEND、IRECV、WAIT 和 collective；
- 保留 run 段的原始消息边界；
- 是 WSE plan 对拍、demand 编译和 BookSim 五档实验的共同输入。

三种 CCDG 的关系：

```text
Raw = setup + run 原始事件
Trim-only = run 原始事件
Compact = run 事件 + burst 折叠
```

## 5. Stage 4：质量闸门

## 5.1 改造前产物

```text
RUN_DIR/quality_gate.txt
```

文本形式：

```text
# 质量闸门 ... 载体: trimonly_16ranks_global.ccdg
loop_time_sec = 0.00184314
barrier_anchored_ranks = 16 / 16
anchor_span_vs_loop_avg_pct = 0.1%
comm_bytes_conserved = true
proc_grid = 4 by 4 by 1 MPI processor grid
dir_classes = 8  dir_off_lattice = 0
```

含义：

- 检查每个 rank 的 run 窗口是否被 BARRIER 正确锚定；
- 检查 compact 前后字节守恒；
- 检查 `K×K×1` 处理器网格；
- short 模式检查消息是否只在二维面邻居/PBC 接缝之间传递。

## 5.2 Phase 0 改造后新增产物

```text
RUN_DIR/wse_plan_validation.json
```

输出形式：

```json
{
  "passed": true,
  "threshold": 0.02,
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
  "kspace_collectives_covered": true,
  "kspace_collectives": [
    {
      "operation": "ALLREDUCE",
      "bytes": 768,
      "plan_count": 16,
      "ccdg_count": 16,
      "passed": true
    }
  ]
}
```

含义：

- 比较 `scope="run"` 的源码 plan 与 trim-only CCDG；
- 点对点消息按二维方向比较消息数和字节数；
- 每个方向允许的字节偏差不超过 2%；
- Kspace collective 必须能在 CCDG 中找到足够数量的同类型、同大小记录。

`quality_gate.txt` 同时新增一行：

```text
wse_plan_vs_ccdg = PASS
(p2p per-direction bytes <= 2%; Kspace collectives covered)
```

因此 Phase 0 后，Stage 4 不只验证 CCDG 自身，还验证源码 emitter 是否漏消息。

## 6. Stage 5：BookSim 五档实验

Phase 0 前后均继续使用 `trimonly.ccdg`，尚未直接读取 `wse_plan.json`。

## 6.1 free

```text
RUN_DIR/
├── booksim_free/trimonly_<R>ranks_global.ccdg
└── booksim_free_result.txt
```

结果文本形式：

```text
booksim_free: ranks=16 mesh=4x4 mode=booksim steps=1 gating=free
cfg=<path>/ccdg_mesh_booksim_free_....cfg
完成: cycles=... sent=... recv=... unresolved=0
```

含义：不做编译和 EST 门控的原始 trace 基线。

## 6.2 free_fold

```text
RUN_DIR/
├── booksim_free_fold/trimonly_<R>ranks_global.ccdg
└── booksim_free_fold_result.txt
```

含义：只做 PBC fold rank 重标号，不删除同步，不使用 EST。

## 6.3 hb / cerebras / ilv

每档产生：

```text
RUN_DIR/
├── booksim_<tag>/
│   ├── demand_<R>ranks_demand.ccdg
│   └── demand_<R>ranks_demand.est
├── demand_<tag>_plan.log
└── booksim_<tag>_result.txt
```

其中 `<tag>` 是：

```text
hb
cerebras
ilv
```

Demand CCDG 仍是 JSON：

```json
{
  "num_ranks": 16,
  "nodes": ["去除运行时同步并按编译计划重建的节点"],
  "cross_rank_edges": ["编译后的依赖"]
}
```

EST 是文本，每行表示节点最早释放周期：

```text
<node_id> <est_cycles>
1024 350
1025 417
```

含义：

- hb：PBC fold 后的 happens-before 正确性档；
- cerebras：增加 XY stage 和串行 wavelet 约束；
- ilv：在 cerebras 基础上使用 interleave wavelet。

BookSim 的详细文件保存在：

```text
booksim2/results/
├── <stem>.cfg
├── <stem>.log
└── <stem>_stats.txt
```

`_stats.txt` 主要形式：

```text
total_sim_cycles = 1168438
compute_cycles = ...
blocked_cycles = ...
congestion_cycles = ...
sched_wait_cycles = ...
```

## 7. Stage 6：最终评估

Phase 0 前后格式不变：

```text
RUN_DIR/
├── evaluation.txt
└── evaluation.json
```

`evaluation.txt` 是人类可读摘要：

```text
档             cycles  cyc/iter  steps/s blocked% cong% sched% idle% unres verdict
free         1,168,438 1,168,438 1,712    ...                         0    PASS
OVERALL: PASS
```

`evaluation.json` 是机器可读结果：

```json
{
  "run_dir": "...",
  "ranks": 16,
  "capture_steps": 1,
  "noc_ghz": 2.0,
  "compute_cap_ops": 25000000000.0,
  "real_loop_time_sec": 0.00184314,
  "carrier": "trimonly",
  "simmode": "free",
  "sims": {
    "free": {
      "total_cycles": 1168438,
      "cycles_per_iter": 1168438,
      "timesteps_per_sec": 1712,
      "packets_sent": 1732,
      "packets_recv": 1732,
      "unresolved": 0,
      "verdict": "PASS"
    }
  },
  "comparisons": []
}
```

PASS 条件：

```text
unresolved == 0
packets_sent == packets_received
所有已运行档位均 PASS
```

## 8. Phase 0 前后的产物差异总结

### 改造前

```text
RUN_DIR/
├── in.lammps
├── lammps.log
├── log.lammps
├── dumpi-*-NNNN.bin
├── dumpi-*.meta
├── trace_<R>ranks_global.ccdg
├── compact_<R>ranks_global.ccdg
├── trimonly_<R>ranks_global.ccdg
├── quality_gate.txt
├── booksim_*/
├── demand_*_plan.log
├── booksim_*_result.txt
├── evaluation.txt
└── evaluation.json
```

### Phase 0 改造后

保留以上文件，并新增：

```text
RUN_DIR/
├── wse_plan.rankNNNN.jsonl
├── wse_plan.kspace.rankNNNN.jsonl
├── wse_plan.json
├── wse_plan_merge.log
└── wse_plan_validation.json
```

新增产物的核心作用：

```text
per-rank JSONL
  = LAMMPS 每个 rank 直接导出的源码通信需求

wse_plan.json
  = 合并后的全局源码通信计划

wse_plan_validation.json
  = 源码 plan 与 DUMPI/CCDG 的逐方向对拍证据
```

最终关系：

```text
DUMPI/CCDG：描述实际发生的 MPI 行为，当前仍服务 BookSim
WSE plan：描述算法层通信需求，作为 Phase 1 编译器输入候选
validation：证明两条支路在现有实验场景下消息数和字节数一致
```
# `run_noc_pipeline.sh` 流程与输入输出

## 1. 脚本用途

该脚本完成一次完整实验：

```text
实验参数
  → 生成 LAMMPS 输入
  → 真实 MPI 运行 LAMMPS
  → DUMPI trace + 源码级 WSE plan
  → 三种 CCDG
  → 数据质量检查
  → 五种 BookSim 对照档位
  → 汇总性能和正确性结果
```

运行入口：

```bash
./run_noc_pipeline.sh <short|long> <ranks> <cu|h2o|lialocl> <atoms> [free|demand|both]
```

示例：

```bash
./run_noc_pipeline.sh short 16 lialocl 2688 both
./run_noc_pipeline.sh long 16 lialocl 2688 free
```

参数含义：

- `short|long`：短程力，或带 PPPM Kspace/FFT 的长程力；
- `ranks`：MPI rank 数，同时也是 BookSim PE 数，必须是完全平方数；
- `cu|h2o|lialocl`：原子体系；
- `atoms`：目标原子数，脚本会选择最接近的合法复制规模；
- `free|demand|both`：选择执行基线、编译档或全部五档。

每次运行创建独立目录：

```text
runs/pipeline/<mode>_<system>_<atoms>a_<ranks>r_<timestamp>/
```

后文简称为 `RUN_DIR`。

## 2. Stage 0：参数和环境检查

### 输入

- 命令行参数；
- `install/bin/lmp`；
- `install/lib/libdumpi.so`；
- `dumpi2ccdg/dumpi2ccdg`；
- `booksim2/ccdg_demand.py`；
- `booksim2/run_ccdg_mesh.sh`；
- WSE plan 开启时的 merge 和 validate 工具。

### 处理

- 检查 mode、体系、原子数和模拟档位是否合法；
- 要求 `ranks=K²`，得到 BookSim 方形 mesh 边长 `K`；
- 创建唯一的 `RUN_DIR`；
- long 模式选择 `coul/long`，short 模式选择 `coul/cut`。

### 输出

- `K`：二维 mesh 的边长；
- `RUN_DIR`：本次实验全部中间文件和结果的根目录。

这一阶段不产生实验数据，只保证后续阶段使用同一个合法配置。

## 3. Stage 1：生成 `in.lammps`

入口是脚本内的 `make_input()`。

### 输入

- `SYSTEM`：Cu、H2O 或 LiAlOCl；
- `MODE`：short 或 long；
- `NATOMS`：目标原子数；
- `RANKS` 和 `K`；
- `CAPTURE_STEPS`：真实 LAMMPS 运行步数，默认 1；
- cases 目录中的 LiAlOCl data 文件或 Cu EAM 势文件。

### 处理

- 选择与目标原子数最接近的 lattice/replicate 规模；
- 固定处理器分解为 `processors K K 1`；
- 设置周期边界、原子类型、质量、电荷和势函数；
- short 使用截断库仑势；
- long 的 H2O/LiAlOCl 使用 `kspace_style pppm 1.0e-4`；
- Cu 使用 EAM，不经过 Kspace，因此 Cu 的 short/long 实际相同；
- 不执行 minimize，只执行 NVE `run CAPTURE_STEPS`。

### 输出

- `RUN_DIR/in.lammps`：LAMMPS 主输入；
- `RUN_DIR/data.LiAlOCl_nvt_charge`：LiAlOCl 数据文件；
- `RUN_DIR/Cu_u3.eam`：Cu 势文件；
- `ACTUAL_ATOMS`：实际生成的原子数。

`in.lammps` 描述的是“计算什么”，后续 LAMMPS、DUMPI 和 WSE plan 都基于同一
份输入，因此不同通信模型之间可以对齐。

## 4. Stage 2：运行 LAMMPS，同时捕获两类通信描述

### 输入

- Stage 1 的 `in.lammps` 和辅助数据文件；
- MPI rank 数；
- DUMPI 动态库；
- `WSE_PLAN_CAPTURE`，默认值为 1。

### 处理

脚本通过 `mpirun` 启动真实 LAMMPS：

```text
LD_PRELOAD=libdumpi.so
DUMPI_OUTDIR=RUN_DIR
LAMMPS_WSE_PLAN=RUN_DIR/wse_plan
```

同一次 LAMMPS 运行产生两条描述通信的支路：

1. DUMPI 在 MPI 层截获实际调用；
2. LAMMPS emitter 在 CommBrick、Grid3d、FFT 和 PPPM 源码位置输出语义 plan。

运行结束后，脚本读取 LAMMPS 的 `Loop time`，并调用
`merge_wse_plan.py` 合并所有 rank 分片。

### 输出

- `lammps.log`：stdout/stderr，包含 Loop time 和错误信息；
- `log.lammps`：LAMMPS 自身日志，包含处理器网格、thermo 等；
- `dumpi-*.bin`：每 rank 的 MPI 事件二进制流；
- `dumpi-*.meta`：DUMPI trace 的 rank/file 索引；
- `wse_plan.rankNNNN.jsonl`：每 rank 的 CommBrick 源码 plan；
- `wse_plan.kspace.rankNNNN.jsonl`：long/PPPM 的 Kspace/FFT plan；
- `wse_plan.json`：合并后的全局源码通信计划；
- `wse_plan_merge.log`：分片数、run 消息数、collective 数和字节数；
- `LOOP_TIME`：真实 LAMMPS run 段耗时。

两类通信描述的含义不同：

- DUMPI 回答“实际发生了哪些 MPI 调用”；
- WSE plan 回答“LAMMPS 哪个算法阶段需要向哪个 rank 传多少数据”。

当前 BookSim 仍使用 CCDG；WSE plan 在 Phase 0 中用于对拍，Phase 1 才会成为
新编译器输入。

## 5. Stage 3：把 DUMPI 转换为三种 CCDG

三次调用 `dumpi2ccdg`，输入都是 Stage 2 的 `dumpi-*.bin/meta`。

### 5.1 Raw CCDG

输出：

```text
trace_<R>ranks_global.ccdg
ccdg_gen.log
```

含义：

- 保留 setup 和 run 的原始 MPI 事件；
- 信息最完整，文件最大；
- 适合追溯问题，不直接作为主 BookSim 载体。

### 5.2 Compact CCDG

设置：

```text
CCDG_TRIM_SETUP=1
CCDG_COMPACT=1
```

输出：

```text
compact_<R>ranks_global.ccdg
ccdg_compact_gen.log
```

含义：

- 删除初始化阶段，只保留 run 窗口；
- 合并连续通信 burst；
- 用于检查压缩前后通信字节是否守恒；
- 节点少，适合观察压缩效果。

### 5.3 Trim-only CCDG

设置：

```text
CCDG_TRIM_SETUP=1
CCDG_COMPACT 未设置
```

输出：

```text
trimonly_<R>ranks_global.ccdg
ccdg_trimonly_gen.log
```

含义：

- 删除 setup，但不折叠 run 内的消息；
- 保留每条通信的边界和依赖关系；
- 是质量检查、demand 编译和 BookSim 五档实验的共同输入。

三者关系：

```text
Raw = setup + run 原始事件
Compact = run 事件 + burst 折叠
Trim-only = run 原始事件
```

## 6. Stage 4：质量闸门

### 输入

- `lammps.log` 和 `log.lammps`；
- compact CCDG 的生成日志；
- trim-only CCDG 及其生成日志；
- `wse_plan.json`。

### 检查内容

1. 每个 rank 是否都用 run 首尾 BARRIER 锚定；
2. BARRIER span 与 LAMMPS Loop time 的偏差；
3. compact 是否保持通信总字节数不变；
4. trim-only 的首节点是否都是 BARRIER；
5. 是否残留 setup BCAST；
6. 每个 rank 的 SEND/ISEND 数量和节点类型分布；
7. WSE plan 与 trim-only CCDG 的逐方向消息数和字节数是否一致；
8. Kspace collective 是否能在 CCDG 中找到对应记录；
9. LAMMPS 处理器网格是否严格为 `K×K×1`；
10. short 模式是否只包含二维面邻居和 PBC 接缝方向。

### 输出

- `quality_gate.txt`：所有检查的简要证据；
- `wse_plan_validation.json`：每个方向的 plan/CCDG 消息数、字节数和偏差；
- 终端 PASS、警告或错误信息。

### 失败语义

- WSE plan 对拍失败、处理器网格错误、short 出现非法方向：硬失败，退出码 2；
- BARRIER 锚定或字节守恒证据不完整：当前打印警告，实验继续；
- span 偏差是诊断信息，不直接决定 BookSim 正确性。

## 7. Stage 5：编译并注入 BookSim

所有档位都来自同一次捕获生成的 `trimonly.ccdg`，因此档位间差异来自放置和调度
机制，而不是不同的 LAMMPS 运行。

### 7.1 `free`

输入：

- 原始 `trimonly.ccdg`；
- PE 计算能力 `CCDG_COMPUTE_CAP`。

处理：

- 不做 rank 重标号；
- 不生成 EST；
- 保留 CCDG 中的同步和依赖，由 PE 就绪后自由注入网络。

输出：

- `booksim_free/trimonly_*.ccdg`；
- `booksim_free_result.txt`；
- BookSim 生成的 cfg、log 和 stats。

含义：原始 trace 载体的反应式基线。

### 7.2 `free_fold`

输入：

- `trimonly.ccdg`。

处理：

- `ccdg_demand.py --fold-only`；
- 只做 PBC fold rank 重标号；
- 不删除同步，不使用 EST。

输出：

- `booksim_free_fold/trimonly_*.ccdg`；
- `booksim_free_fold_result.txt`。

含义：单独观察 PBC 放置变化的影响。

### 7.3 `hb`

输入：

- `trimonly.ccdg`；
- `--fold-pbc`；
- 计算能力参数。

处理：

- 提炼通信需求；
- 应用 PBC fold；
- 生成保证 happens-before 正确性的无同步 demand CCDG；
- 生成节点最早释放时间 EST。

输出：

- `booksim_hb/demand_*_demand.ccdg`；
- `booksim_hb/demand_*_demand.est`；
- `demand_hb_plan.log`；
- `booksim_hb_result.txt`。

含义：需求驱动编译的正确性基准。

### 7.4 `cerebras`

在 hb 参数上增加：

```text
--cerebras
```

输出结构与 hb 相同，目录和日志标签为 `cerebras`。

含义：增加 Cerebras 风格 XY stage 和 wavelet 调度约束，用于测量阶段串行化及
barrier 代价。

### 7.5 `ilv`

在 cerebras 参数上增加：

```text
--wavelet-mode interleave
```

输出结构与 hb 相同，目录和日志标签为 `ilv`。

含义：相位交错版本，用于与串行 wavelet 调度比较。

### 模式选择

- `free`：只运行 free；
- `demand`：运行 hb、cerebras、ilv；
- `both`：运行 free、free_fold、hb、cerebras、ilv。

脚本前部仍保留了一组较早的 `inject_free()`/`inject_demand()` 定义；后面的五档
版本重新定义了同名函数，Shell 实际执行后面的版本。理解当前实验时应以后面的
五档实现为准。

## 8. Stage 6：结果汇总

### 输入

- 每个实际运行档位的 `booksim_<tag>_result.txt`；
- 结果中引用的 BookSim cfg、stats 和 log；
- rank 数、捕获步数、Loop time 和计算能力。

### 提取指标

- `total_sim_cycles`：BookSim 完成整个 CCDG 的总周期；
- `cycles_per_iter`：总周期除以捕获步数；
- `timesteps_per_sec`：按 2 GHz NoC 换算的吞吐；
- `compute_cycles`：PE 执行计算节点的累计周期；
- `blocked_cycles`：等待通信依赖的累计周期；
- `congestion_cycles`：因网络/注入队列拥塞等待的周期；
- `sched_wait_cycles`：等待 EST 释放的周期；
- `idle`：总 PE-cycle 账本中未被以上项目占用的比例；
- `packets_sent/received`：网络收发包数；
- `unresolved`：仿真结束后仍未解决的跨 rank 依赖。

账本关系：

```text
ranks × makespan
  = compute + blocked + congestion + sched_wait + done/idle
```

### PASS 条件

- `unresolved == 0`；
- 如统计到 packets，则 `packets_sent == packets_received`；
- 所有已运行档位都 PASS，整体才是 PASS。

### 输出

- `evaluation.json`：机器可读的完整指标和档位对比；
- `evaluation.txt`：人类可读的摘要；
- 脚本退出码：
  - `0`：全部档位通过；
  - `1`：至少一个 BookSim 档位失败；
  - `2`：输入、捕获、转换或质量闸门错误。

档位间自动比较：

- free → free_fold：PBC fold 放置影响；
- free_fold → hb：需求驱动编译净收益；
- hb → cerebras：XY stage/wavelet 串行化代价；
- cerebras → ilv：串行和交错调度差异。

## 9. 完整产物树

```text
RUN_DIR/
├── in.lammps
├── lammps.log
├── log.lammps
├── dumpi-*.bin
├── dumpi-*.meta
├── wse_plan.rankNNNN.jsonl
├── wse_plan.kspace.rankNNNN.jsonl
├── wse_plan.json
├── wse_plan_merge.log
├── wse_plan_validation.json
├── trace_<R>ranks_global.ccdg
├── compact_<R>ranks_global.ccdg
├── trimonly_<R>ranks_global.ccdg
├── ccdg_gen.log
├── ccdg_compact_gen.log
├── ccdg_trimonly_gen.log
├── quality_gate.txt
├── booksim_free/
├── booksim_free_fold/
├── booksim_hb/
├── booksim_cerebras/
├── booksim_ilv/
├── booksim_<tag>_result.txt
├── demand_<tag>_plan.log
├── evaluation.txt
└── evaluation.json
```

并非每次运行都会产生全部 BookSim 目录；具体取决于 `free|demand|both`。
BookSim runner 生成的 `.cfg`、`.log` 和 `_stats.txt` 实际保存在
`booksim2/results/`，`booksim_<tag>_result.txt` 中记录了对应 cfg 路径；
Stage 6 根据该路径找到并解析这些文件。

## 10. 阅读结果时的核心关系

```text
in.lammps
  └→ 决定原子体系、处理器网格和通信需求

DUMPI trace
  └→ dumpi2ccdg
       ├→ raw：追溯完整 MPI
       ├→ compact：验证折叠和字节守恒
       └→ trim-only：五档共同实验载体

WSE plan
  └→ 与 trim-only 对拍，证明源码 emitter 没有漏消息

trim-only
  ├→ free / free_fold
  └→ demand compiler → hb / cerebras / ilv

BookSim logs/stats
  └→ evaluation.txt/json
```

因此，定位问题时应按以下顺序检查：

1. `lammps.log`：LAMMPS 是否正确运行；
2. `ccdg_*_gen.log`：裁剪和转换是否正确；
3. `quality_gate.txt`：输入数据能否用于实验；
4. `wse_plan_validation.json`：源码 plan 是否与 trace 对齐；
5. `demand_<tag>_plan.log`：编译器是否生成合法计划；
6. `booksim_<tag>_result.txt` 和 stats：网络仿真是否完成；
7. `evaluation.txt`：最终性能和正确性结论。
