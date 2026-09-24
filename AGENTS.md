# LAMMPS MPI Trace → NoC Simulation Pipeline

## 项目概述

本项目构建了一个**闭环仿真流水线**，用于将 LAMMPS 的 MPI 通信 trace 转换为通信-计算依赖图（CCDG），并驱动 BookSim 2.0 进行片上网络（NoC）仿真。

### 流水线架构

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│  LAMMPS      │     │  DUMPI       │     │  dumpi2ccdg  │     │  BookSim 2.0 │
│  MPI 仿真    │ ──→ │  Trace 捕获  │ ──→ │  CCDG 转换   │ ──→ │  NoC 仿真    │
└──────────────┘     └──────────────┘     └──────────────┘     └──────────────┘
     │                     │                     │                     │
     │ 运行 LAMMPS         │ LD_PRELOAD 注入      │ libundumpi 解析     │ CCDGTrafficManager
     │ 生成 MPI 事件       │ 捕获通信 trace       │ 构建依赖图          │ 驱动网络仿真
     │                     │                     │                     │
     ↓                     ↓                     ↓                     ↓
   in.lammps        dumpi-*.bin/.meta       trace_*.ccdg         results.csv
```

---

## 目录结构

```
/work1/jiangtao/lammps_trace/
│
├─【构建与安装】──────────────────────────────────────────────────────
├── lammps-src/                # LAMMPS 源码（29 Aug 2024，KSPACE 相位插桩、remap_3d 转置 marker 在此）
├── lammps-build/              # 原生 MPI 构建（OpenMPI，KSPACE+MANYBODY）→ make install 到 install/
├── lammps-build-smpi/         # SMPI 版构建（smpirun 专用，产物 lmp，重链接后须 strip DF_1_PIE）
├── install/                   # 原生版安装目录
│   ├── bin/lmp                # 原生 MPI 版 LAMMPS（DUMPI 流水线用）
│   ├── bin/dumpi2ascii        # DUMPI trace 转文本
│   └── lib/libdumpi.so        # DUMPI 捕获库（LD_PRELOAD 注入）
├── sst-dumpi/                 # libdumpi 源码（DUMPI 捕获/回放库）
│
├─【转换工具】────────────────────────────────────────────────────────
├── dumpi2ccdg/                # DUMPI trace → CCDG 转换器（dumpi2ccdg.cpp + 可执行）
├── booksim2/                  # BookSim 2.0 NoC 仿真器（CCDGTrafficManager 驱动）
│
├─【算例】cases/               # LAMMPS 输入算例
│   ├── lialocl_coul/          # LiAlOCl 基座（84 原子，data 文件被其他变体复用）
│   ├── lialocl_coul_{1344,2688}/   # replicate 变体（SMPI 强扩展主算例，-v STEPS）
│   ├── lialocl_coul_2688_nokspace/ # 无 Kspace 变体（coul/cut）
│   ├── lialocl_coul_{4,8,16,32}r/  # 各 rank 数旧变体
│   ├── h2o_coul/              # 2 型 LJ+Coul 合成体系（fcc）
│   ├── eam_cu_1k/             # Cu/EAM 势（Cu_u3.eam）
│   ├── lj_bench/, lj_coul_*   # LJ / LJ+Coul 基准
│   └── smpi_test/             # SMPI 冒烟用例
│
├─【SMPI mesh 强扩展实验】───────────────────────────────────────────
├── run_kspace_profile.sh      # Kspace 操作级拆分（单例/batch，产物 runs/<ts>_batch/）
├── simgrid_traces/            # 2D mesh 平台 XML（platform_mesh_{2x2,2x4,4x4,4x8}.xml）
│   ├── gen_mesh_platform.py   #   平台生成（1.42Gf/6.8GBps/10ns）
│   └── strip_pie_flag.py      #   lmp 重链接后清 DF_1_PIE（SMPI 必需）
├── runs/
│   ├── smpi_mesh_{2688,1344}/ # 主实验：LAMMPS 原生 breakdown（真实时钟）
│   ├── smpi_fft_split/        # FFT 计算/转置通信分离（analyze_fft_split.py）
│   ├── smpi_trace_kspace/     # Paje+host-speed 相位解剖（attribute_kspace.py）
│   └── <时间戳>_batch/        # run_kspace_profile.sh 的 batch 产物
│
├─【DUMPI→CCDG→SimGrid 验证流水线】──────────────────────────────────
├── run_short_lmp.sh           # 无 Kspace 版（入参: <rank数> <cu|h2o|lialocl> <原子数> [步数]）
├── run_long_lmp.sh            # 含 Kspace(PPPM) 版（同参，Cu/EAM 两者等价）
├── run_lmp_ccdg.sh            #   共用引擎：生成输入→DUMPI 捕获→CCDG→验证
├── ccdg_simgrid_validate.py   #   SimGrid DAG 验证驱动（默认自动校准频率，SIMGRID_CPU_FREQ 为固定口径）
├── validate_simgrid.py        # CCDG SimGrid DAG 验证（全 4/8/16/32 ranks 批量版）
├── runs/pipeline/             # 流水线产物 <mode>_<体系>_<原子数>a_<rank数>r_<时间戳>/
│   └── in.lammps / lammps.log / log.lammps / dumpi-*.bin|meta /
│       trace_<N>ranks_global.ccdg / simgrid_dag/ / validation_result.txt
├── run_trace_capture.sh       # 手动 DUMPI trace 捕获（老接口）
│
├─【分析与验证脚本】──────────────────────────────────────────────────
├── ccdg2dag.py / ccdg2simgrid.py / compute_tsimgrid.py / run_dag_simulation.py  # CCDG→SimGrid DAG
├── ccdg_critical_path.py      # 关键路径分析
├── validate_iterative.py / validate_kspace_autocal.py / validate_lialocl_weakscaling.py
│
├─【文档】────────────────────────────────────────────────────────────
├── AGENTS.md                  # 本文件（项目总览与快速上手）
├── Kspace_FFT实验手册.md       # FFT 分离实验完整方案/结果/坑清单（上下文交接文档）
├── Kspace_profile.md          # Kspace 操作级拆分方法
├── dumpi2ccdg_guide.md        # CCDG 转换详细指南
├── CCDG实验报告.md / CCDG_vrf_tasklist.md / tasklist.md
├── .qoder/skills/smpi-mesh-experiment/  # 项目 Skill（SMPI mesh 实验规范环境，自动触发）
├── validation_report_simgrid.json       # CCDG 验证报告
└── traffic_pattern/           # 通信模式分析脚本与结果
```

---

## 环境依赖

### 已安装组件

| 组件 | 路径 | 说明 |
|------|------|------|
| LAMMPS | `install/bin/lmp` | MPI 分子动力学仿真器 |
| LAMMPS (SMPI) | `lammps-build-smpi/lmp` | SMPI 重编译版（链接 libsimgrid），用于 SimGrid 虚拟平台重放 |
| libdumpi | `install/lib/libdumpi.so` | DUMPI 动态捕获库 |
| dumpi2ascii | `install/bin/dumpi2ascii` | Trace 转文本工具 |
| dumpi2ccdg | `dumpi2ccdg/dumpi2ccdg` | Trace → CCDG 转换器 |
| BookSim | `booksim2/booksim` | NoC 网络仿真器 |
| SimGrid 3.30 | 系统包（smpirun/smpicc） | SMPI 仿真平台，驱动 2D mesh 虚拟平台实验 |

### 编译命令

```bash
# 编译 dumpi2ccdg
cd dumpi2ccdg
make

# 编译 BookSim
cd booksim2/src
make
```

---

## 实验流程

### 阶段一：Trace 捕获

**脚本**: [run_trace_capture.sh](file:///work1/jiangtao/lammps_trace/run_trace_capture.sh)

```bash
# 捕获 4-rank trace
sudo ./run_trace_capture.sh 4 /work1/jiangtao/lammps_trace/cases/lj_bench

# 捕获 16-rank trace
sudo ./run_trace_capture.sh 16 /work1/jiangtao/lammps_trace/cases/lj_bench
```

**输出**: `runs/trace_{N}ranks_{timestamp}/dumpi-*.bin/.meta`

### 阶段二：CCDG 转换

**工具**: [dumpi2ccdg.cpp](file:///work1/jiangtao/lammps_trace/dumpi2ccdg/dumpi2ccdg.cpp)

```bash
# 转换 4-rank trace
cd dumpi2ccdg
./dumpi2ccdg ../runs/trace_4ranks_20260721_201258 > ../runs/trace_4ranks_20260721_201258/trace_4ranks_global.ccdg

# 转换 16-rank trace
./dumpi2ccdg ../runs/trace_16ranks_20260721_201406 > ../runs/trace_16ranks_20260721_201406/trace_16ranks_global.ccdg
```

**输出**: `trace_{N}ranks_global.ccdg`

**详细说明**: 参见 [dumpi2ccdg_guide.md](file:///work1/jiangtao/lammps_trace/dumpi2ccdg_guide.md)

### 阶段三：单配置仿真

**配置文件**: [ccdg_lammps_4x4.cfg](file:///work1/jiangtao/lammps_trace/booksim2/ccdg_lammps_4x4.cfg)

```bash
cd booksim2
./booksim ccdg_lammps_4x4.cfg
```

**输出**: `lammps_4ranks_stats.txt`

### 阶段四：多配置批量仿真

**脚本**: [run_configs.sh](file:///work1/jiangtao/lammps_trace/booksim2/run_configs.sh)

```bash
cd booksim2
bash run_configs.sh
```

**输出**: `results/all_results.csv`

---

## 核心设计决策

### 1. Rank → NoC Node 映射

**规则**: rank 号直接作为 NoC node 号（一一映射）

```
Rank 0 → Node 0, Rank 1 → Node 1, ..., Rank N-1 → Node N-1
```

**约束**: `_num_ranks == _nodes`（CCDG 的 rank 数必须等于 NoC 节点数）

### 2. 频率比例转换

**公式**: `freq_ratio = cpu_freq_ghz / noc_freq_ghz`

**应用**: 将 CCDG 中的 CPU 计算周期转换为 NoC 周期

```cpp
pe.remaining_cycles = node.compute_cycles / _freq_ratio;
```

### 3. 依赖解析机制

**msg_id 映射**: 每个数据包分配唯一 ID，映射到跨 rank 边

```cpp
// 发送时记录映射
_msg_id_to_edge_idx[msg_id].push_back(edge_idx);

// 接收时解析依赖
_cross_edges[edge_idx].resolved = true;
_advancePE(dst_rank);  // 唤醒阻塞的 PE
```

### 4. PE 状态机

| 状态 | 行为 |
|------|------|
| `PE_COMPUTE` | 递减 `remaining_cycles`，完成后前进 |
| `PE_BLOCKED` | 等待跨 rank 依赖解决 |
| `PE_DONE` | 所有节点处理完毕 |

---

## LAMMPS 单次迭代通信清单（强扩展分析基准）

适用场景：lj/cut/coul/long + PPPM + `fix nve` + brick 分解（本项目 lialocl_coul 系列用例）。源码依据：`lammps-src/src/verlet.cpp`（迭代主循环）、`src/comm_brick.cpp`（原子通信）、`src/KSPACE/pppm.cpp`（kspace 通信）、`src/grid3d.cpp`（网格交换）。

### 迭代时间线与通信位置

```
┌─ COMM ──── exchange() ①  →  borders() ②
├─ NEIGH ─── neighbor->build()          （本地计算，无通信）
├─ PAIR ──── pair->compute()            （本地计算，无通信）
├─ BOND ──── 成键相互作用（lialocl 用例无，跳过）
├─ KSPACE ── PPPM::compute()            ← 通信最密集
│             ③ reverse_comm（ρ 网格幽灵）
│             ④ brick2fft + poisson（FFT 转置 ×2 + 3D FFT）
│             ⑤ forward_comm（E 场网格幽灵）
│             ⑥ fieldforce（本地）
│             ⑦ MPI_Allreduce（能量/virial）
├─ COMM ──── reverse_comm() ⑧           （newton on，力归约）
├─ MODIFY ── fix nve 积分                （本地）
└─ OUTPUT ── thermo 输出（每 thermo N 步一次，含 Allreduce）
```

### 通信清单

| # | 通信 | 相位 | 类型 | 通信模式 | 局部/全局 |
|:--:|---|---|---|---|---|
| ① | `CommBrick::exchange()` 原子迁移 | Comm | P2P（Irecv+Send） | 6 面邻居 | **局部** |
| ② | `CommBrick::borders()` 幽灵原子 | Comm | P2P + 尺寸协商 Sendrecv/Allreduce | 6 面邻居 | **局部** |
| ③ | PPPM `reverse_comm` 电荷密度 ρ | Kspace | P2P swap | 6 面邻居 | **局部** |
| ④ | PPPM FFT 转置（Remap/FFT3d） | Kspace | 平面交换，全排列对等通信 | 与所有相关 rank | **全局** |
| ⑤ | PPPM `forward_comm` 电场 E | Kspace | P2P swap | 6 面邻居 | **局部** |
| ⑦ | PPPM `MPI_Allreduce` 能量+virial | Kspace | 集合通信 | world 全员 | **全局** |
| ⑧ | `CommBrick::reverse_comm()` 力 | Comm | P2P swap | 6 面邻居 | **局部** |
| — | thermo 热力学量归约 | Output | `MPI_Allreduce` | world 全员 | **全局** |

### 关键机制特征

1. **所有 P2P 通信只发生在 brick 分解的 6 面邻居之间**（①②③⑤⑧），消息大小随 rank 数增加而缩小，最终进入延迟主导区；**唯一的真·全局通信集中在 kspace**（④⑦）与稀疏的 thermo 输出。
2. **通信域不可区分**：Grid3d/FFT3d/Remap 全部使用同一个 `world` 通信域，DUMPI 的 comm 字段无法区分 kspace 内的局部/全局通信，必须依赖相位追踪标记。
3. **FFT 转置（④）是强扩展瓶颈**：数据需在 brick 分解与 FFT 分解间重排，通信量无法随 rank 数摊薄。实测（2688 原子，4→32r）：全局通信 wall 时间 0.55→5.52 s（×10，超线性），占总通信 49.8%→58.8%。
4. **网格交换多轮退化**：当某维度 brick 网格面数跌破 PPPM stencil ghost 宽度时，`Grid3d::setup_comm_brick` 需多轮 swap 接力数据，所有 rank 陪跑相同轮数（含零字节空消息）。实测 16r 时 GRID_REVERSE 事件数/rank/步从 ~12 增至 ~18，其中约 31% 为零字节消息。
5. **邻居表频率影响通信结构**：`neigh_modify every 1` 时①②每步发生；生产运行用 `every 10 + delay` 会稀疏化①②，kspace 通信占比进一步升高。
6. **消息尺寸阈值法在强扩展下失效**：所有消息随 rank 数变小，2KB 阈值误分类率高，必须用源码级相位标记分类（见下节）。

### 相位追踪（Kspace 全局通信精确捕捉）

已实现源码级标记机制，用代码位置而非消息特征归类通信：

- **`src/timer.h/.cpp`**：rank 0 将相位切换/子相位标记写入 `phase_trace.csv`（CLOCK_MONOTONIC ns，与 DUMPI wall time 同基准），通过环境变量 `LAMMPS_PHASE_TRACE=<路径>` 启用；相位行为 END 标记，`KSPACE_*` 行为子阶段 START 标记
- **`src/KSPACE/pppm.cpp`**：`compute()` 中 5 个标记——局部：`KSPACE_GRID_REVERSE`/`KSPACE_GRID_FORWARD`；全局：`KSPACE_FFT`/`KSPACE_REDUCE`；另有 `KSPACE_FORCE`
- **`dumpi2ccdg.cpp`**：每个通信节点携带 `wall_time_sec`/`wall_duration_sec`（统计信息已改走 stderr，stdout 纯 JSON）
- **`traffic_pattern/analyze_phase.py`**：按 wall time 将 CCDG 通信节点归位到相位窗口与 kspace 子阶段；单实验模式或 `--batch` 多实验强扩展趋势对比
- **注意**：输入脚本需 `timer normal`（默认）以上级别，`timer loop/off` 会静默丢失相位标记

```bash
# 批量强扩展趋势分析（每个 run_dir 需含 *.ccdg + phase_trace.csv）
python3 traffic_pattern/analyze_phase.py --batch runs/trace_4ranks_* runs/trace_8ranks_* \
    runs/trace_16ranks_* runs/trace_32ranks_* \
    --csv traffic_pattern/results/trend.csv
```

**统计口径提醒**：`analyze_phase.py` 输出的是"通信内部结构占比"（全 rank MPI 调用时长求和）；LAMMPS log 的 `MPI task timing breakdown %total` 是"相位 wall 时间/总 wall 时间"（含计算）。两者口径不同，报告中不可混用；强扩展下 log 的 Kspace %total 单调上升（计算缩小、通信不缩），而通信内部结构占比受非 kspace 通信"地板"影响可能非单调。

---

## 配置参数说明

### CCDGTrafficManager 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `sim_type` | - | 必须设置为 `ccdg` |
| `ccdg_file` | - | CCDG 文件路径 |
| `cpu_frequency_ghz` | 2.5 | CPU 频率 |
| `noc_frequency_ghz` | 1.0 | NoC 频率 |
| `flit_size_bytes` | 8 | 每个 flit 的字节数 |

### BookSim 拓扑参数

| 参数 | 说明 |
|------|------|
| `topology` | 网络拓扑（mesh/torus/fattree） |
| `k` | 维度基数（mesh: 每行/列节点数） |
| `n` | 维度数（mesh: 2=二维） |
| `num_vcs` | 虚拟通道数 |
| `vc_buf_size` | 每个 VC 的 buffer 大小 |

### 常用配置组合

| Rank 数 | 拓扑配置 | 节点数 |
|---------|----------|--------|
| 4 | `topology=mesh, k=2, n=2` | 4 |
| 16 | `topology=mesh, k=4, n=2` | 16 |
| 8 | `topology=mesh, k=2, n=3` | 8 |

---

## 验证方法

### CCDG 验证

**脚本**: [validate_ccdg.sh](file:///work1/jiangtao/lammps_trace/validate_ccdg.sh)

```bash
./validate_ccdg.sh runs/trace_4ranks_20260721_201258/trace_4ranks_global.ccdg
```

**检查项**:
- JSON 格式正确性
- 节点 ID 唯一性
- 边引用有效性
- Send-Recv 匹配
- 统计信息合理性

### 仿真验证

**检查项**:
- 仿真是否正常结束（非死循环）
- 发送/接收数据包数相等
- 延迟和吞吐率在合理范围

---

## 常见问题

### Q1: rank 数与 NoC 节点数不匹配

**现象**: 数据包发送到错误位置或无法注入

**解决**: 确保配置文件中的拓扑节点数等于 CCDG 的 rank 数

### Q2: 仿真死锁或未完成

**可能原因**:
- 1x4 mesh 等退化拓扑（k=1 导致所有节点映射到同一 router）
- 依赖解析不完整（msg_id 未找到）

**解决**: 使用 2x2、4x4 等正方形拓扑

### Q3: CCDG 生成失败

**可能原因**:
- DUMPI trace 文件损坏或不完整
- MPI 事件匹配失败

**解决**: 使用 `dumpi2ascii` 检查原始 trace

### Q4: 频率比例设置

**建议**:
- CPU 频率通常为 2.5-3.0 GHz
- NoC 频率根据仿真需求调整（0.001-1.0 GHz）
- 频率比越大，计算节点执行时间越短

---

## 参考文档

| 文档 | 路径 | 说明 |
|------|------|------|
| CCDG 转换指南 | [dumpi2ccdg_guide.md](file:///work1/jiangtao/lammps_trace/dumpi2ccdg_guide.md) | 详细转换流程和数据结构 |
| 任务清单 | [tasklist.md](file:///work1/jiangtao/lammps_trace/tasklist.md) | 完整任务列表和结果分析 |
| CCDG 验证脚本 | [validate_ccdg.sh](file:///work1/jiangtao/lammps_trace/validate_ccdg.sh) | CCDG 正确性验证工具 |

---

## 快速上手命令

```bash
# 1. 捕获 trace（4-rank）
sudo ./run_trace_capture.sh 4 /work1/jiangtao/lammps_trace/cases/lj_bench

# 2. 转换为 CCDG
cd dumpi2ccdg
./dumpi2ccdg ../runs/trace_4ranks_* > ../runs/trace_4ranks_*/trace_4ranks_global.ccdg

# 3. 验证 CCDG
cd ..
./validate_ccdg.sh runs/trace_4ranks_*/trace_4ranks_global.ccdg

# 4. 运行仿真
cd booksim2
./booksim ccdg_lammps_4x4.cfg

# 5. 查看结果
cat lammps_4ranks_stats.txt
```

---

## CCDG 验证（SimGrid DAG）

### 验证目的

使用 SimGrid 4.1 DAG 仿真验证 `dumpi2ccdg` 提取的通信-计算依赖图（CCDG）的正确性。通过对比仿真时间 $T_{\text{simgrid}}$ 与 LAMMPS 真实执行时间 $T_{\text{real}}$ 来判定，验证阈值为 $E_{\text{total}} \le 5\%$。

### 验证流水线

```
DUMPI Trace ──► dumpi2ccdg ──► CCDG JSON ──► SimGrid DAG C++ ──► 编译运行 ──► T_sim vs T_real
```

**脚本**: [validate_simgrid.py](file:///work1/jiangtao/lammps_trace/validate_simgrid.py)

```bash
# 使用环境变量指定 CPU 频率运行（推荐用 2.8 GHz 一致频率）
SIMGRID_CPU_FREQ=2.80 python3 validate_simgrid.py

# 或使用自动校准（从 rank=4 校准频率用于所有 ranks）
python3 validate_simgrid.py

# 仅验证指定 rank 数
python3 validate_simgrid.py --ranks 4,8,16,32
```

### 关键改进

| 改进项 | 说明 |
|--------|------|
| **ALLREDUCE 递归加倍** | 用 log₂(N) 步顺序通信取代原始 N-1 个并发 P2P 发送，rank=32 Comm 数降 84% |
| **BCAST 二分树** | root 按 1→2→4→8→... 步进发送，减少网络争用 |
| **一致 CPU 频率** | 使用 2.8 GHz（Intel Xeon Gold 5418Y 全核睿频）而非逐 rank 自动校准，避免通信开销被吸收进计算速度 |

### 验证结果（CPU 频率 = 2.8 GHz 一致）

| Ranks | 活动数 | Comm 数 | $T_{\text{real}}$ | $T_{\text{simgrid}}$ | $E_{\text{total}}$ | 状态 |
|:-----:|:------:|:-------:|:-----------------:|:--------------------:|:------------------:|:----:|
| 4 | 23,940 | 2,681 | 0.7085 s | 0.6898 s | **2.64%** | PASS |
| 8 | 71,548 | 9,088 | 0.3447 s | 0.3332 s | **3.35%** | PASS |
| 16 | 148,740 | 22,795 | 0.2155 s | 0.2138 s | **0.78%** | PASS |
| 32 | 308,579 | 53,753 | 0.1482 s | 0.1492 s | **0.69%** | PASS |

**结论**: dumpi2ccdg 在 rank=4, 8, 16, 32 上全部验证通过。

### 核心发现

1. **CPU 频率逐 rank 校准的陷阱**：公式 `freq = max_cycles / T_real` 隐含假设所有时间用于计算，但 `T_real` 包含 MPI busy-wait 开销。Rank 数越多，通信占比越大，校准频率越虚高（rank=32 达 3.92 GHz）。使用 CPU 真实物理频率（2.8 GHz）更为合理。

2. **ALLREDUCE 建模精度影响最大**：该操作是 LAMMPS 中最频繁的集合通信（占通信节点 16-17%），原始 N-1 模型在 rank=32 时产生 213,280 个通信活动，远超实际的递归加倍算法（34,400 个），导致 25.68% 的误差。

### 详细文档

| 文档 | 路径 |
|------|------|
| 完整实验报告 | [CCDG实验报告.md](file:///work1/jiangtao/lammps_trace/CCDG实验报告.md) |
| 验证结果 JSON | [validation_report_simgrid.json](file:///work1/jiangtao/lammps_trace/validation_report_simgrid.json) |
| 验证脚本 | [validate_simgrid.py](file:///work1/jiangtao/lammps_trace/validate_simgrid.py) |
| 任务列表 | [CCDG_vrf_tasklist.md](file:///work1/jiangtao/lammps_trace/CCDG_vrf_tasklist.md) |

---

## SMPI 重放（2D mesh 虚拟平台）：Kspace 占比与 Rank 数关系

### 方法

不用 DAG 仿真，直接让真实 LAMMPS 二进制跑在 SimGrid 的 2D mesh 虚拟平台上（SMPI replay）：SMPI 版 `lmp`（`lammps-build-smpi/lmp`，用 smpicxx 重编译、链接 libsimgrid）由 `smpirun` 通过 dlopen 加载，MPI 调用由 SMPI 按 mesh 网络模型计时。**LAMMPS 自带 Timer 基于 `MPI_Wtime`，在 SMPI 下返回仿真时钟**，因此 log 中的 `MPI task timing breakdown` 与 `Loop time` 直接就是 mesh 上的仿真结果。

**关键运行条件**（缺一即失败）：

| 条件 | 原因 |
|------|------|
| `SMPI_PRIVATIZATION=0` | privatization（dlopen 副本）会使 `MPI_COMM_WORLD` 绑定为 NULL，报 `MPI_ERR_COMM` |
| 二进制清除 `DF_1_PIE` 标志 | glibc 拒绝 dlopen 带 PIE 标志的可执行文件；用 `simgrid_traces/strip_pie_flag.py` 清除（每次重新链接后需重做） |
| `-Wl,--export-dynamic-symbol=main` | SMPI dlopen 后需从动态符号表找 `main` |

注意：log 中有**两个** `MPI task timing breakdown`，第一个属于 minimize 阶段，第二个才是 MD run——提取数据时勿取错。

**已知限制**：`LD_PRELOAD` 拦截对 dlopen 加载的程序无效（符号直接绑定 libsimgrid）。逐调用计时 wrapper 需静态链接进二进制（`simgrid_traces/smpi_walltime_static.c`，定义 `MPI_*` 调 `PMPI_*`，产出 `lmp_wt`），当前 SMPI fork 模型下 rank 解析未调通，仅能依赖 LAMMPS 相位级 timer。

### 平台与命令

平台由 [gen_mesh_platform.py](file:///work1/jiangtao/lammps_trace/simgrid_traces/gen_mesh_platform.py) 生成（`platform_mesh_{2x2,2x4,4x4,4x8}.xml`），当前参数：**host 1.42 Gflops，链路 6.8 GB/s，延迟 10 ns**；rank r → host_r 坐标 `(r % KX, r // KX)`，BFS 最短路显式路由。

```bash
cd cases/lialocl_coul_2688   # 或 lialocl_coul_1344
SMPI_PRIVATIZATION=0 smpirun \
    -platform /work1/jiangtao/lammps_trace/simgrid_traces/platform_mesh_4x4.xml \
    -np 16 /work1/jiangtao/lammps_trace/lammps-build-smpi/lmp \
    -in in.lammps -v STEPS 10
```

### 实验结果（Kspace %total 随 rank 数变化）

算例：LiAlOCl（lj/cut/coul/long + PPPM + nve），10 步 MD run。

**2688 原子**（FFT 网格 40×24×36），日志在 `runs/smpi_mesh_2688/`：

| ranks | mesh | Loop time (s) | Kspace %total | Kspace avg (s/rank) | Comm %total | Pair %total |
|:-----:|:----:|:-------------:|:-------------:|:-------------------:|:-----------:|:-----------:|
| 4  | 2×2 | 0.1495 | 72.5% | 0.1083 | 18.7% | 6.2% |
| 8  | 2×4 | 0.1596 | 78.3% | 0.1250 | 17.5% | 2.3% |
| 16 | 4×4 | 0.2913 | 80.7% | 0.2350 | 16.9% | 0.7% |
| 32 | 4×8 | 1.5252 | **92.9%** | 1.4172 | 6.3% | 0.1% |

**1344 原子**（replicate 4 2 2，FFT 网格 36×20×20），日志在 `runs/smpi_mesh_1344/`：

| ranks | mesh | Loop time (s) | Kspace %total | Kspace avg (s/rank) | Comm %total | Pair %total |
|:-----:|:----:|:-------------:|:-------------:|:-------------------:|:-----------:|:-----------:|
| 4  | 2×2 | 0.0551 | 71.8% | 0.0396 | 18.3% | 6.3% |
| 8  | 2×4 | 0.0910 | 76.6% | 0.0697 | 18.7% | 2.1% |
| 16 | 4×4 | 0.2252 | 82.0% | 0.1847 | 15.6% | 0.5% |
| 32 | 4×8 | 0.8053 | **88.7%** | 0.7143 | 9.9% | 0.1% |

### 核心结论（第一阶段目标）

1. **Kspace 占比随 rank 数单调上升**：2688 原子 72.5% → 92.9%，1344 原子 71.8% → 88.7%；两种规模趋势一致，是结构性规律而非算例巧合。
2. **Loop time 超线性恶化且增量几乎全来自 Kspace**：2688 原子 4→32 ranks Loop time ×10.2，Kspace 绝对耗时 ×13；1344 原子更差（×14.6）——系统越小越先撞通信墙。
3. **瓶颈是带宽与多跳争用而非延迟**：延迟从 1 μs 降到 10 ns 未扭转趋势；带宽减半（12.5 → 6.8 GB/s）使 32r Loop time ×2.47（4r 仅 ×1.41）——FFT 转置的全对全流量在多跳链路上汇聚，是全局通信超线性恶化的直接证据。
4. **局部通信（Comm 相位）扩展性显著优于 kspace 内全局通信**：Comm 绝对耗时 4→32 仅 ×3.4，占比从 ~19% 掉到 6-10%，被 Kspace 的爆炸彻底淹没；1344 在 32r 时 Comm 占比反升，系小 brick 触发网格幽灵多轮交换（空消息陪跑）所致。
5. **计算能力影响甚微**：host speed 减半（2.85 → 1.42 Gf）后 Pair 绝对耗时翻倍但占比 <6%，强扩展区间内系统为通信受限（communication-bound）。

### 平台参数敏感性对照（2688 原子）

| 平台配置 | 4r Loop/Kspace% | 32r Loop/Kspace% |
|----------|:---------------:|:----------------:|
| 2.85 Gf / 12.5 GB/s / 1 μs | 0.106 s / 72.5% | 0.617 s / 87.7% |
| 1.42 Gf / 6.8 GB/s / 10 ns | 0.149 s / 72.5% | 1.525 s / 92.9% |

### 待办（第二阶段）

- 修复 `lmp_wt` 静态 wrapper 的 rank 解析（SMPI fork 模型下 `MPI_Comm_rank` 只在 fork 前的父进程调用），实现逐 MPI 调用仿真计时
- 结合 `phase_trace.csv` 与 `analyze_phase.py`，把 Kspace 内部再拆为全局/局部通信的精确占比
