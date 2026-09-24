---
name: lmp-booksim
description: >-
  LAMMPS MPI trace → compact CCDG → BookSim 2D mesh NoC 仿真的完整操作规范。
  当用户要求跑 LAMMPS trace 捕获/CCDG 转换、把 CCDG 注入 BookSim mesh、
  做 free/WSE/demand 三档注入对比、SimGrid 单步验证、或分析 BookSim cycles
  结果（unresolved、blocked、cycles/iter、算力扫描）时使用本 skill。
  也适用于用户提到 run_short_lmp.sh / run_long_lmp.sh / run_ccdg_mesh.sh /
  ccdg_demand.py / compact_v2 / trimonly 等"trace 到 NoC"链路工作。
---

# lmp-booksim：LAMMPS trace → CCDG → BookSim mesh 注入

## 当前 WSE 编译入口（Phase 3）

纯短程论文口径优先使用静态 WSE 编译链，不再从 DUMPI trace 反推 demand：

```bash
# 单点：PLAN/COST -> 静态波前程序 -> BookSim 验收
python3 wse_compiler.py PLAN COST booksim2/ccdg_lammps_4x4.cfg \
  --wse-fast-profile --compute-capability 2.5e10 -o OUT
booksim2/run_wse_program.sh OUT.program.json

# 12 点标准矩阵：16/64/256 × fast/IQ × fold on/off
python3 run_wse_phase3.py
```

- 设计、口径与结果：`WSE_compiler_design.md`
- 结果：`experiments/wse_phase3/summary.json`、`matrix.csv`、`per_stage.csv`、
  `link_utilization.csv` 和 `heatmap_*.svg`
- fold 默认开启，`--no-fold` 关闭；实验必须显式指定 `--compute-capability`
- 12 点必须同时满足 compiler/BookSim 偏差 ≤1%、branch 全送达和 unresolved=0

## 顶层入口

一条命令完成"真实 LAMMPS 跑一次 → DUMPI 捕获 → CCDG 生成 → 质量闸门 → 四档注入 → evaluation"：

```bash
./run_noc_pipeline.sh <short|long> <rank数> <cu|h2o|lialocl> <原子数> [模式=both]
# 模式: free（trace 自由注入上界）| demand（三档 demand 编译）| both（默认四档同源对比）
# 四档: free / hb（happens-before 正确性保证）/ cerebras（XY stage+wavelet）/ ilv（相位交错对照）
```

- **无 SimGrid**：数据质量闸门 = BARRIER 锚定（run 首尾 barrier 对 vs Loop time）+ comm_bytes 守恒 + 每 rank 首节点 BARRIER；注入后硬校验 unresolved==0 且 sent==recv
- **Z 压缩载体**：`processors K K 1`（Z 邻居本地化，Cerebras 范式），闸门校验处理器网格与方向集（只允许 2D 面邻居 ±1 / PBC 接缝 ±(k-1)）
- 产物：`runs/pipeline/<mode>_<体系>_<原子数>a_<rank数>r_<时间戳>/`，评估输出 `evaluation.txt/.json`
- 退出码：0=PASS，2=流水线/闸门错误
- 旧入口 `run_short_lmp.sh`/`run_long_lmp.sh`（含 SimGrid）保留但不再维护

## 历史五档口径（trace-demand 基线归档）

| 档 | 语义 | 时钟语义 |
|---|---|---|
| free | trace 自由注入（logical rank 空间），同步节点保留 | 反应式上界 |
| free_fold | free + paper III-E PBC 交叉排布重标号 | 隔离 fold 对反应式档的影响 |
| hb | demand 编译 + happens-before 约束，**fold-PBC 载体**（消费点 est ≥ feeder 到达） | 正确性保证基准 |
| cerebras | hb + X/Y 阶段重排 + 全局 stage barrier + b+1 相位串行轮转 | 论文机制忠实复现 |
| ilv | hb + stage barrier + 相位交错格点（无串行轮转） | 隔离相位轮转代价 |

五档用于历史基线和 `编译式 WSE - trace-driven demand` 净收益对照，不再作为当前
WSE 主实现。只有算例、算力和捕获来源一致时才可作强结论。

## 关键产物命名（同一 run 目录内）

| 文件 | 含义 |
|---|---|
| `trace_<N>ranks_global.ccdg` | 原始 CCDG（含 setup） |
| `compact_<N>ranks_global.ccdg` / `compact_v2_*.ccdg` | 裁剪+折叠，仅迭代段；**v2 = BARRIER 锚定窗口版（当前标准口径）** |
| `trimonly_<N>ranks_global.ccdg` | 裁剪不折叠，demand 编译器输入载体 |
| `unrolled_<k>steps.ccdg` | BookSim 多步注入用展开版 |
| `simgrid_dag/` + `validation_result.txt` | SimGrid DAG 单步验证 |

## 手动注入 BookSim

```bash
cd booksim2
./run_ccdg_mesh.sh <ccdg文件> [超时秒=3600] [步数=1] [wse_gating=0] [phase_width=64] [strip_width=2]
```

- 要求 `num_ranks` 为完全平方数（4/9/16/25/64/256…），mesh 取 k=√N, n=2
- 频率口径：NoC 2.0 GHz，折算墙钟 = cycles × 0.5 ns
- 算力：`CCDG_COMPUTE_CAP`（ops/s，默认 2.5e10）或 `CCDG_COMPUTE_RATE`（ops/cycle，**必须带小数点**，整数 0 走错 Assign 分支）
- 注入反压：`CCDG_INJECT_QDEPTH`（flit 数，0=auto）
- WSE 门控：第 4-6 参或 `wse_gating=1`
- CSV 追加写入 `booksim2/results/ccdg_mesh_results.csv`（含 unresolved、cycles_per_iter、blocked 等列）

## 结果判定（每次跑完必查）

1. **unresolved == 0 且 sent == recv**：否则注入有问题，结果不可用
2. `cycles_per_iter` = total_cycles / 步数；`timesteps_per_sec = 1e9/(cycles_per_iter×0.5ns)`
3. 与真实 Loop time 交叉核对：Σcompute/rank ÷ 校准频率 ≈ Loop time（100–104% 区间为正常）
4. SimGrid 精度验证：
   ```bash
   python3 ccdg_simgrid_validate.py <run_dir> <num_ranks> <compact_v2_*.ccdg>
   ```
   自动校准含 `SIMGRID_FREQ_CAP_GHZ=4.0` 物理上限检查；固定频率用 `SIMGRID_CPU_FREQ`。阈值 E_total ≤ 5%。

## 口径陷阱（易踩坑）

- **旧 compact（Loop-time 墙钟窗口版）会裁掉 halo P2P**，64r 时 0 SEND/rank，"完美线性"结论失效；必须用 `compact_v2`（BARRIER 锚定）
- **跨载体对比无意义**：同为 256r，不同捕获批次（包数不同）的 demand vs free 结论方向可以相反；只做同源对比
- **demand 旧口径的 blocked=0 是删除语义不是满足语义**：不带 HB 约束时 68%（16r）/36%（256r）依赖存在"消费早于到达"，makespan 是乐观下界；论文级结论必须用 hb 档
- **bypass 规划序 ≠ 发射链序**：cerebras 模式 pop-scan 后必须按 plan_t 重排链再 emit，否则 est 不可达（16r 实测差 83%），已内置修复
- **PBC 接缝**：K proc 一维时周期缝邻居相距 k-1 跳（16r b=3、256r b=15），是 torus→mesh 映射的真实代价，不是 Z 投影残留；闸门只允许面邻居 ±1 与接缝 ±(k-1)
- **SimGrid 256r 退化**：强扩展下 compute 占比趋零时自动校准频率无物理意义（曾得 0.005 GHz），这是移除 SimGrid 闸门的原因之一
- log 里有**两个** `MPI task timing breakdown`，第一个是 minimize，第二个才是 MD run

## 历史数据与已确立结论

具体数字、旧实验复核、对比基线：读 `references/baseline_results.md`。当前确立的口径（Z 压缩 + fold-PBC 载体，2026-08-31）：

- **五档同源结果（全部 unresolved=0、sent==recv、断言全绿）**：16r free 175,046 / free_fold 174,119 / hb 163,844 / cerebras 165,892 / ilv 165,903；256r free 151,285 / free_fold 146,286 / hb 137,867 / cerebras 162,694 / ilv 162,729
- **需求驱动净收益（正确性口径）**：+5.9%（16r）/ +5.8%（256r，fold 载体），稳定
- **fold-PBC（paper III-E 交叉排布）的反向发现**：free 档 fold 赚（256r +3.3%，接缝 b=15 长飞行消失）；hb 档 fold 中性（±0.6%）；**cerebras 档 fold 净亏（256r hb→cerebras 从 +14.1% 恶化到 +18.0%）**——单播投影下 fold 把"少数接缝长流"换成"全体主流 b=1→b=2 链路占用翻倍"，被 stage 串行放大。论文 "PBC is free" 依赖多播原语（b=2 波前 fan-out 覆盖 3 tiles ≠ 2 倍单播成本）；复现该结论需要 MCAST 支持（见待办）
- **Cerebras stage 串行代价**：16r +1.3~1.6%，256r +14~18%——代价随 k 增长，与 fold 与否弱相关；主因是串行化放大链路占用而非接缝本身
- **相位串行 vs 交错几乎无差**（≤0.08% 全场景一致）：代价在 barrier 不在轮转
- **est↔实测偏差**：全部档 ≤0.02%（链序一致修复后）
- compact_v2 单步 SimGrid 验证 4/4 ≤1%（历史口径，SimGrid 已从主流水线移除）
- long/short cycles 比 ~1.6–3.9×，Kspace 全局通信是强扩展瓶颈

## 需求驱动编译细节

需要单独操作或调参数时读 `booksim2/ccdg_demand_workflow.md`。编译器 `ccdg_demand.py` 档位：默认 = hb（正确性保证）；`--no-hb` = 旧乐观下界；`--cerebras` = XY stage 串行化 + barrier + wavelet；`--stage-barrier=global|rank|none`、`--wavelet-mode=serial|interleave`。BookSim 二进制与 `run_ccdg_mesh.sh` 零改动，复用 `ccdg_schedule_file` 通道。
