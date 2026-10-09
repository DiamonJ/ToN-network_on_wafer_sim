# LAMMPS trace → CCDG → BookSim NoC 仿真：复现包

本包包含完整的源码、脚本、文档、算例与一个样例 run，用于复现
"LAMMPS MPI 通信 trace → 通信-计算依赖图 (CCDG) → BookSim 2D mesh NoC 注入 →
需求驱动 (Cerebras 范式) 编译对比" 的整条实验流水线。

参考论文（未随包分发，请自行获取）：
- Santos et al., "Breaking the Molecular Dynamics Timescale Barrier Using a
  Wafer-Scale System", SC 2024, IEEE/ACM.（wavelet / b+1 相位 / fold-PBC / 多播机制出处）
- BookSim 2.0: Jiang et al., Stanford.（随包源码，含本项目扩展）

当前纯短程主入口是编译式 WSE：

```bash
python3 run_wse_phase3.py
```

它复现 16/64/256 ranks × WSE-fast/IQ × fold on/off 的 12 点矩阵。设计和实测结论见
`WSE_compiler_design.md`，结果见 `experiments/wse_phase3/`。下文 DUMPI→CCDG 五档
流水线作为历史 trace-demand 基线保留。

## 1. 目录结构

```
├── README_REPRODUCE.md      本文件
├── smoke_test.sh            一键构建 + 冒烟验证
├── AGENTS.md                项目总览（架构、口径、坑清单；必读）
├── run_noc_pipeline.sh      ★ 主入口：LAMMPS 捕获 → CCDG → 五档注入 → evaluation
├── wse_compiler.py          ★ 当前短程 WSE 静态波前编译器
├── run_wse_phase3.py        ★ 12 点 WSE 评估、per-stage 指标和链路热图
├── WSE_compiler_design.md   当前架构、复现协议、结果与边界
├── experiments/wse_phase3/  Phase-3 机器可读结果和热图
├── run_short_lmp.sh / run_long_lmp.sh / run_lmp_ccdg.sh   旧流水线（SimGrid 验证口径，保留）
├── ccdg_unroll.py           单步 CCDG → 多步展开（BookSim 多步注入用）
├── validate_*.py, ccdg_simgrid_validate.py 等              SimGrid DAG 验证（可选）
├── dumpi2ccdg/              DUMPI trace → CCDG 转换器（C++ 源码）
├── booksim2/                BookSim 2.0（含 CCDGTrafficManager、WSE 扩展）
│   ├── src/                 C++ 源码（ccdg_trafficmanager.cpp 为核心扩展）
│   ├── ccdg_demand.py       ★ 需求驱动编译器（hb/cerebras/fold-PBC 各档）
│   ├── ccdg_planner.py / ccdg_scheduler.py                时槽表与调度工具
│   ├── run_ccdg_mesh.sh     CCDG → BookSim mesh 注入
│   ├── ccdg_demand_workflow.md  需求驱动编译工作流与结论（必读 §11–§12）
│   └── results/             历史实验 CSV/stats（结果证据）
├── sst-dumpi/               DUMPI 捕获/回放库源码（libdumpi / libundumpi）
├── lammps-src/              LAMMPS 29Aug2024-Update1 源码
│                            （已插桩：timer.cpp/h 的 LAMMPS_PHASE_TRACE、
│                              KSPACE/pppm.cpp 相位标记）
├── cases/                   LAMMPS 输入算例（lialocl / cu / h2o / lj 等）
├── simgrid_traces/          2D mesh 平台 XML 生成（SMPI 实验，可选）
├── traffic_pattern/         通信模式/相位分析脚本
├── docs/skill/              lmp-booksim skill（操作规范 + 历史基线数据）
└── runs/pipeline/short_lialocl_2688a_16r_20260831_195028/
                             ★ 样例 run（16 rank：DUMPI trace + 各档 CCDG +
                               BookSim 结果 + evaluation），可不跑 LAMMPS 直接复用
```

## 2. 依赖

| 组件 | 版本要求 | 用途 |
|---|---|---|
| g++ | ≥ 9（支持 C++20，SimGrid DAG 用） | 全部 C++ 构建 |
| OpenMPI | ≥ 4.x（`mpirun --allow-run-as-root` 可用） | LAMMPS 捕获运行 |
| python3 | ≥ 3.8，标准库即可 | 编译器/评估脚本 |
| GNU make / autoconf | 系统自带 | sst-dumpi、booksim 构建 |
| SimGrid 3.30+（可选） | `apt install simgrid` 或源码 | 仅旧流水线 SimGrid 验证需要 |
| 硬件 | 多核 x86_64 Linux；16 rank 捕获需 ≥16 逻辑核（`--oversubscribe` 已加） | |

## 3. 构建（~15 分钟）

假设解压后的根目录为 `$PKG`（脚本使用相对自身定位的绝对路径，无需安装到系统）：

```bash
export PKG=$(pwd)

# ① sst-dumpi: libdumpi(LD_PRELOAD 捕获) + libundumpi(读取) → 安装到 $PKG/install
cd $PKG/sst-dumpi
./configure --prefix=$PKG/install --enable-libdumpi CC=mpicc CXX=mpic++
make -j$(nproc) && make install

# ② dumpi2ccdg
cd $PKG/dumpi2ccdg && make            # DUMPI_PREFIX 默认 ../install

# ③ BookSim（含 CCDGTrafficManager）
cd $PKG/booksim2/src && make -j$(nproc)

# ④ LAMMPS（原生 MPI 版，已含相位插桩；KSPACE+MANYBODY 足够本项目算例）
#    注意: WITH_PNG/JPEG/FFMPEG/GZIP 必须显式关闭（与原始构建一致），
#    否则系统 libpng 误探测会导致链接失败 png_set_longjmp_fn
cd $PKG/lammps-src
mkdir build && cd build
cmake ../cmake -D PKG_KSPACE=yes -D PKG_MANYBODY=yes \
      -D BUILD_MPI=yes -D WITH_PNG=no -D WITH_JPEG=no -D WITH_FFMPEG=no -D WITH_GZIP=no \
      -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=$PKG/install
make -j$(nproc) && make install       # → $PKG/install/bin/lmp
```

## 4. 一键冒烟验证

```bash
cd $PKG && bash smoke_test.sh
```

依次执行：三件套构建检查 → 用样例 run 的 CCDG 直接注入 BookSim（无需 MPI）→
跑完整 `run_noc_pipeline.sh short 16 lialocl 2688 both`（真实 LAMMPS 捕获 + 五档编译/注入 +
评估）。两条都 PASS 即环境正确。

## 5. 历史五档基线流水线用法

```bash
cd $PKG
./run_noc_pipeline.sh <short|long> <rank数> <cu|h2o|lialocl> <原子数> [模式=both]
# rank 数须为完全平方数（4/9/16/25/64/256…），mesh 取 k=√N
# 示例（本项目主算例）:
./run_noc_pipeline.sh short 16  lialocl 2688 both
./run_noc_pipeline.sh short 256 lialocl 2688 both     # 约 10 分钟
```

五档注入（同一次捕获内对比，跨 run 目录对比无意义）：

| 档 | 语义 |
|---|---|
| free | trace 自由注入（logical rank 空间），反应式上界 |
| free_fold | free + 论文 III-E PBC 交叉排布重标号 |
| hb | demand 编译 + happens-before 正确性约束（默认口径） |
| cerebras | hb + X/Y 阶段串行化 + stage barrier + b+1 相位轮转 |
| ilv | hb + barrier + 相位交错（隔离轮转代价的对照档） |

产物在 `runs/pipeline/<mode>_<体系>_<原子数>a_<rank数>r_<时间戳>/`：
`evaluation.txt / evaluation.json`（评估）、`quality_gate.txt`（质量闸门证据）、
`trace|compact|trimonly_*ranks_global.ccdg`、`booksim_*/`（各档）、
`wse_plan.json`（Plan 发射）、`static_cost_estimate.json`（T1/C1 估算）、
`wse_phase1.program.json / .replay.ccdg / .replay.est / .report.json / .acceptance.json`
（WSE 编译产物，仅 short 模式默认生成）。

以下子步骤中，**前三项已内置在流水线中**，后两项为独立脚本：

**① WSE Plan 发射**（流水线内置，`WSE_PLAN_CAPTURE=1` 默认启用）

```bash
# 流水线自动执行：mpirun 捕获时 LAMMPS_WSE_PLAN 已设 → merge → validate
# 如需独立运行：
LAMMPS_WSE_PLAN=./wse_plan mpirun -np 16 lmp -in in.lammps   # 每 rank 写 wse_plan.rank<N>.jsonl
python3 merge_wse_plan.py ./wse_plan --expected-ranks 16 -o wse_plan.json
python3 validate_wse_plan.py wse_plan.json trimonly_16ranks_global.ccdg --threshold 0.02
```

**② 静态 T1/C1 估算 + PMU 验证**（估算始终执行，PMU 默认关闭需 `COMPUTE_CAPTURE=1`）

```bash
# 静态估算：流水线自动执行
python3 estimate_lammps_cost.py in.lammps -o static_cost_estimate.json
# PMU 捕捉 + 对拍：流水线中通过环境变量开启
COMPUTE_CAPTURE=1 ./run_noc_pipeline.sh short 16 lialocl 2688 both
# 或独立运行：
bash profile_lammps_flops.sh 16 in.lammps 100 lmp              # 100 = 测量步数（扣除 run 0 基线）
python3 validate_static_cost.py static_cost_estimate.json wse_plan.json compute_profile.json
```

**③ WSE 编译**（流水线内置，`WSE_COMPILE=1` 默认启用，仅 short 模式）

```bash
# 流水线自动执行：compile → BookSim 验收
# 如需独立运行（需先生成 wse_plan.json 和 static_cost_estimate.json）：
python3 wse_compiler.py wse_plan.json static_cost_estimate.json ccdg_lammps_4x4.cfg \
  --wse-fast-profile --compute-capability 2.5e10 -o wse_phase1
bash booksim2/run_wse_program.sh wse_phase1.program.json       # BookSim 验收（cycles≤1%，counts 守恒）
```

**④ 多原子规模批量实验**（独立脚本，不在流水线内）

```bash
python3 run_wse_phase3.py --all
# 产物：experiments/wse_phase3/ → per_stage.csv, link_utilization.csv, 热力图 PNG/SVG
```

**⑤ 注入率扫图**（独立脚本，不在流水线内）

```bash
cd booksim2 && make -C src && python3 sweep_queue_time.py       # → sweep_results/inj_vs_queue_time.png
```

## 6. 验收标准与参考数字

每次 run 必须满足（evaluation.txt 中 OVERALL: PASS 即自动包含）：

1. 质量闸门：BARRIER 锚定 = ranks/ranks、comm_bytes 守恒、方向集 ⊆ {面邻居 ±1, PBC 接缝 ±(k−1)}
2. 每档 BookSim：`unresolved == 0` 且 `packets_sent == packets_recv`
3. 编译器断言：happens-before 违例 = 0、stage_inversions = 0（ cerebras/ilv 档）、
   fold 后 max_b ≤ 2
4. est ↔ 实测 makespan 偏差 ≤ 1%

参考数字（本包样例 run，2.5 GHz host 口径 cap=2.5e10 ops/s；**绝对值随机器浮点性能浮动，
以相对关系与断言为准**）：

| 载体 | 16r (4×4) | 256r (16×16) |
|---|---:|---:|
| free | 175,046 | 151,285 |
| free_fold | 174,119 | 146,286 |
| hb | 163,844 | 137,867 |
| cerebras | 165,892 | 162,694 |
| ilv | 165,903 | 162,729 |

已确立结论（详见 `docs/skill/references/baseline_results.md` 与
`booksim2/ccdg_demand_workflow.md` §11–§12）：需求驱动净收益（正确性口径）+5.8~5.9%；
Cerebras stage 串行代价 16r ~+1.3%、256r +14~18%；单播投影下 fold-PBC 对 cerebras 档
净亏（依赖多播原语，见 §12）；相位串行 vs 交错 ≤0.08%。

## 7. 已知坑（务必读）

- `dumpi2ccdg` 的环境变量开关**只查存在性不查值**：要 trimonly 必须 `env -u CCDG_COMPACT`；
  置 `CCDG_COMPACT=0` 仍会触发折叠（流水线已处理，手工转换时注意）。
- `CCDG_COMPUTE_RATE` 必须带小数点（整数 0 走错 Assign 分支）。
- LAMMPS log 有两个 `MPI task timing breakdown`，第一个是 minimize，第二个才是 MD run。
- 旧 `compact_*`（Loop-time 墙钟窗口）会裁掉 halo P2P；当前流水线用 BARRIER 锚定版
  （`ccdg_trimonly_gen.log` 中 "step anchored on BARRIER pair" 即证据）。
- 手工 BookSim 注入时，同一目录放多个 CCDG 会因 TAG=父目录名互相覆盖 stats。
- 256r 强扩展下 compute 占比趋零，SimGrid 自动校准频率会退化到无物理意义值
  （0.005 GHz 量级）——这是 SimGrid 闸门从主流水线移除的原因。

## 8. 环境说明与许可

- 本包为内部研究复现用途整合：BookSim 2.0 遵循其原始许可证（Stanford，见
  `booksim2/src` 内版权声明）；LAMMPS 为 GPL（`lammps-src/LICENSE`）；sst-dumpi 遵循其
  上游许可证（SST/mirror 项目）；其余脚本与文档随项目分发。
- `runs/` 历史数据未随包携带（1.7 GB），样例 run 之外请用主流水线重新生成。
