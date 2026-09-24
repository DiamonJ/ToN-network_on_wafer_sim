# WSE 编译器设计、实验与复现

## 1. 范围

v1 面向 SC'24 论文覆盖的纯短程 MD 稳态迭代。LAMMPS 在 MPI 运行时输出
`wse_plan.<rank>.json`，合并器形成逻辑通信计划；`wse_compiler.py` 将该计划编译为
静态 compute block、forward exchange/borders 多播波前和 reverse reduction 波前，
最后由 BookSim `WSETrafficManager` 回放。

long/PPPM 暂不纳入主结论。kspace 网格 swap 可以沿用波前，但 FFT remap 是全排列通信，
论文没有对应机制，必须作为边界实验单列。

## 2. 数据通路

```text
LAMMPS (N ranks)
  -> wse_plan.<rank>.json
  -> merge_wse_plan.py
  -> wse_plan.json + static_cost_estimate.json
  -> wse_compiler.py
       program.json       语义程序、placement、stage、wavefront、link footprint
       est                worker compute 时间表
       replay.ccdg/est    逐链路分支的 BookSim 可执行回放
       report.json        守恒、断言、周期模型
  -> run_wse_program.sh
       booksim.stats      manager 实测与逐 stage 完成周期
       acceptance.json    1% makespan 闸门
```

编译入口：

```bash
python3 wse_compiler.py PLAN COST booksim2/ccdg_lammps_4x4.cfg \
  --wse-fast-profile --compute-capability 2.5e10 -o OUT
booksim2/run_wse_program.sh OUT.program.json
```

`--no-fold` 禁用 fold-PBC；不传时默认开启。`--wse-fast-profile` 使用当前可执行的
4-byte flit / 3-cycle hop WSE 参数；不传时使用基础 IQ 配置（8-byte flit /
4-cycle hop）。`--compute-capability` 覆盖 cfg 中算力，实验必须显式传值，避免 WSE
与历史基线使用不同算力。

## 3. 编译语义

1. C0/C1 静态 compute block 在每个 worker 上执行，并以全局 compute barrier 收口。
2. 每个方向的逻辑消息按轴、forward/reverse 和 fold 后的物理坐标归组。
3. forward 用根到叶的多播树；reverse 用相反方向的归约树。
4. 每棵树按 `b+1` 轮形成 wavefront，`b` 是该 stage 的最大树深。
5. replay 将树边投影成相邻 tile 的 branch packet，并对 injection/ejection/link
   资源做确定性区间调度。
6. manager 只按编译 release cycle 发射，不从运行时 trace 反推 demand。

核心正确性条件：

- logical message bytes 输入守恒；
- placement 为双射，fold 开启时 PBC 最大跳数受界；
- wavefront 内链路互斥；
- reduction 叶先于根；
- replay command/branch 全部送达；
- compiler makespan 与 BookSim makespan 偏差不超过 1%。

## 4. Phase 3 实验协议

一键入口：

```bash
python3 run_wse_phase3.py
```

矩阵固定为 `{16,64,256} × short × {WSE-fast, IQ} × {fold on,off}`，共 12 点。
三档输入都是 LiAlOCl、单个短程稳态 timestep；64-rank 输入由
`cases/wse_short_64/in.lammps` 真实运行产生。NoC 为 2 GHz，compute capability 固定
为 `2.5e10 ops/s`。结果在 `experiments/wse_phase3/`：

- `summary.json`：完整机器可读结果；
- `matrix.csv`：cycles/iter、steps/s、拥塞、波前及字节；
- `per_stage.csv`：forward/reverse × H/V 的编译周期和 manager 实测完成周期；
- `link_utilization.csv` 与 12 张 `heatmap_*.svg`：方向链路负载；
- `baseline_comparison.csv`：编译式 WSE 对历史 trace-demand 基线的净收益。

## 5. Phase 3 结果

所有 12 点均 PASS：分支送达率 100%，unresolved 为 0，compiler/BookSim 相对误差
为 0.0065%–0.0178%，没有观察到运行时拥塞。WSE-fast fold-on 的
cycles/iter 为 107,553 / 75,194 / 82,937（16/64/256 ranks），对应
18.6k / 26.6k / 24.1k steps/s。

IQ fold-on 为 82,811 / 47,131 / 44,996 cycles，比 WSE-fast 少 23.0% /
37.3% / 45.8%。这不是“IQ 路由器更快”的结论：当前两个 profile 同时改变了
flit width，IQ 的 8-byte flit 抵消了其更长 hop pipeline；它是一项硬件参数敏感性
扫描。若要隔离路由器实现，必须增加同 flit width 的 profile。

fold-on 相对 fold-off 在 WSE-fast 下分别为 -0.02% / -0.33% / +0.81%。
256-rank IQ 下 fold 收益为 +2.37%。因此在编译式多播/归约模型里，PBC 成本已经接近
free（绝对影响不超过 2.4%），修复了历史单播 demand 投影里 fold 反向亏损的诊断；
但不能把该范围外推到 FFT remap。

WSE-fast fold-on 的编译期通信分解（H-forward / V-forward / H-reverse /
V-reverse）：

- 16r：21.3% / 28.7% / 21.3% / 28.7%；
- 64r：13.6% / 36.4% / 13.6% / 36.4%；
- 256r：9.6% / 40.4% / 9.6% / 40.4%。

V 轴 payload 更大，规模增大后成为主导项。热点链路的聚合 flit-cycle /
通信窗口从 16r 的 0.253 降至 64r 的 0.205、256r 的 0.141；完整空间分布见 SVG
热图。`per_stage.csv` 的 manager span 允许不同 stage 尾部重叠，不能直接相加；
需要可加占比时使用 compiler start/completion。

## 6. 历史基线与证据边界

16-rank 同一短程算例、同一 `2.5e10 ops/s` 口径有完整五档历史结果：
free / free_fold / hb / cerebras / ilv。WSE-fast fold-on 相对正确性基准 hb
减少 29.18% cycles，相对 cerebras 减少 29.79%；这就是“编译式 WSE vs
trace-driven demand”的主要净收益。

256-rank 当前只自动匹配到历史 free 基线，净收益为 47.99%；它来自历史 CSV，
不是同一次捕获，证据等级低于 16-rank 五档，只用于趋势参考。64-rank 没有
`2.5e10` 的五档同源历史记录，因此不伪造跨口径比较。旧五档实现和既有结论已归档到
`booksim2/ccdg_demand_workflow.md`，不再作为 WSE 主入口。

## 7. 风险与后续

- fast/IQ 参数扫描耦合了 flit width 与 pipeline，后续应拆成正交敏感性实验。
- 当前公开资料没有 WSE 核精确算力；`--compute-capability` 必须保留并做扫描，
  不能把单点绝对 steps/s 当硬件预测。
- v2 路由器只能作为可选优化；任何接入都要重新通过 credit、死锁和 12 点矩阵。
- long/PPPM 仅可作为机制边界实验，不进入论文短程性能结论。
# Phase-1 WSE Compiler

`wse_compiler.py` 将 LAMMPS 源码级 `wse_plan.json` 编译成确定性的 short-range
WSE program。MPI trace 不参与编译，只保留为输入计划的校准 oracle。

## 输入和输出

```bash
python3 wse_compiler.py \
  runs/.../wse_plan.json \
  runs/.../static_cost_estimate.json \
  booksim2/ccdg_lammps_4x4.cfg \
  -o runs/.../wse_phase1
```

输出：

- `wse_phase1.program.json`：worker placement、compute blocks、H/V stages、
  `(b+1)` phase wavefront、multicast/reduction flow、逐跳 footprint；
- `wse_phase1.est`：`wavefront_id earliest_start_cycle`；
- `wse_phase1.report.json`：守恒、冲突、placement、fold 和 EST 断言。
- `wse_phase1.replay.ccdg/.replay.est`：manager-level multicast 的 BookSim
  回放载体；它是 WSE program 的后端 lowering，不是 MPI trace。

`run_noc_pipeline.sh` 在 short 模式下默认执行该步骤；设置
`WSE_COMPILE=0` 可关闭。编译器直接读取
`booksim2/ccdg_lammps_4x4.cfg`。其中 `k` 是 pipeline 模板值，实际 mesh
尺寸由 `wse_plan.procgrid` 决定；flit、NoC 频率、计算能力和 router pipeline
delay 均使用配置文件值。Pipeline 额外启用 `--wse-fast-profile`：4-byte flit
和 IQ router 可执行的最浅流水（routing=0、VC=1、switch=1，hop stride=3）。

## Phase-1 语义

1. 要求 `Px × Py × 1` 方形 worker grid，完成 Z 压缩。
2. 偶数边长使用论文 III-E fold-PBC 交叉排布，周期邻居不超过 2 hops。
3. CommBrick `dimension=0/1` 分别编译为 H/V stage。
4. 每个 stage 的 `b=max(round)+1`，条带 phase 数为 `b+1`。
5. forward 类 phase 编译为 multicast forest；reverse 类 phase 标记为
   `reduction_op=sum`，并导出 2:1 分支候选节点。
6. `TreeLinkTable` 对方向链路和 destination ejection 进行 flit 级预留；
   H/V stage 和 phase tick 使用 completion barrier 串联。

当前 source plan 只给出每次聚合消息的总字节数，没有 payload lineage，无法证明
两条重叠路径携带同一份数据。因此 footprint 对同一 `(link,depth)` 只预留一个
共享 slot，但 slot 负载取各 payload 之和，保证消息/字节守恒，不会为了得到更低
数字而错误地丢掉 payload。在 LAMMPS emitter 增加 payload/root lineage 后，
slot 负载可无歧义地改为真正 multicast 去重。

`booksim2/run_wse_program.sh` 使用 `sim_type=wse` 和
`WSETrafficManager` 回放 compiler lowering，并生成
`wse_phase1.acceptance.json`。相对误差超过 1% 时命令返回失败。

## 边界

- 仅接受 `forward/reverse/pair_forward/pair_reverse`；FFT、Grid3d 和 collective
  会显式报错。
- 静态 C1 当前作为通信前的全局 compute barrier；细粒度计算/通信重叠需要新的
  source-level compute phase emitter。
- 当前是 traffic-manager 级 v1：compiler 将共享 `(link,depth)` footprint slot
  展开为一次一跳 branch packet，并为每个 wavefront 生成一个单 flit command；
  BookSim 真实执行 DOR、router pipeline、VC、credit 和 ejection。
  路由器内动态复制 flit 的 v2 仍为可选后续项。

## 回归

```bash
python3 -m unittest tests.test_wse_compiler
```

已用真实 16-rank 与 256-rank LiAlOCl short 运行产物完成全断言和 BookSim
端到端验收：16r 为 628853/628860 cycles（误差 0.0011%，192/192 branches），
256r 为 143547/143554 cycles（误差 0.0049%，3840/3840 branches）；
两者 congestion ratio 均为 0。
