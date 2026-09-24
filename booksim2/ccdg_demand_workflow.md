# [历史归档] ccdg_demand.py 需求驱动五档工作流

> 本文保留 free/free_fold/hb/cerebras/ilv 五档的历史实现与实验口径，用作
> trace-driven demand 基线。当前纯短程 WSE 主入口是 `wse_compiler.py` 和
> `run_wse_phase3.py`；设计及 12 点结果见 `../WSE_compiler_design.md`。
> 新实验不要从本文的无 HB 旧口径推导正确性结论。

## 1. 概述

**目标**：按 Cerebras SC24 范式，从 LAMMPS 的**通信需求**（谁给谁发什么、几何方向模式）出发，由编译器在 2D mesh 上排布注入时刻，消除 MPI 运行时同步带来的 blocked cycles。

**核心思想**：

- MPI trace 里的 WAIT/WAITANY/BARRIER 是运行时实现产物，不是物理需求；
- 从 trace 提炼**消息需求实例集**（丢弃时序与同步结构），消息顺序全局自由；
- 时槽表（LinkTable）保证链路无冲突、数据按编排时刻确定性到达；
- 旧实现删除接收等待后得到 `blocked = 0 by construction`；这只是历史乐观下界，
  正确性对照必须使用后续加入 happens-before 的 hb 档；
- BookSim 二进制与 `run_ccdg_mesh.sh` **零改动**（复用 `ccdg_schedule_file` 的 sched_est 门控通道）。

**适用范围**：P2P halo 短程（方向聚类 ≤8 类）。FFT 转置/沿线多播为后续阶段。

## 2. 工作流总览

**已整合进流水线**：`run_short_lmp.sh` / `run_long_lmp.sh` 的 `sim=demand` 模式一条命令完成下列全部阶段（真实 LAMMPS 只跑一次）；所有中间产物都在本次运行目录 `runs/pipeline/<mode>_<sys>_<atoms>a_<r>r_<ts>/` 下，不再用 /tmp。

```
① trace 捕获            run_lmp_ccdg.sh（LAMMPS 真实 1 步 DUMPI trace）
        │
② trimonly CCDG         dumpi2ccdg + CCDG_TRIM_SETUP=1（env -u CCDG_COMPACT）
        │  含 COMPUTE/SEND/ALLREDUCE + WAIT*/RECV/BARRIER 同步节点
        ▼
③ 需求提炼              ccdg_demand.py: extract_demand()
        │  丢弃同步节点，仅留 COMPUTE/SEND/collective 需求链
        │  校验 comm_bytes 守恒 + 方向聚类直方图
        ▼
④ 时槽排布              ccdg_demand.py: plan_flow()
        │  PipelineLinkTable 无冲突时槽表，pe_t/port_t 双时间线
        │  collective burst 整体平移搜索
        ▼
⑤ 无同步 CCDG 生成      ccdg_demand.py: emit()
        │  <prefix>_demand.ccdg（无 cross_rank_edges、无同步节点）
        │  <prefix>_demand.est（"node_id 释放时刻"表）
        ▼
⑥ BookSim 两档注入      booksim_free/（trimonly free 基线）+ booksim_demand/（est 门控）
        │
⑦ 统计与账本            stats 文件：total/blocked/congestion/sched_wait
                        账本守恒：ranks×makespan = compute+blocked+congestion+sched_wait+PE_DONE_idle
```

## 3. 阶段 1：一条命令整合流程（sim=demand）

```bash
cd /work1/jiangtao/lammps_trace
./run_short_lmp.sh 16 lialocl 2688 1 demand   # long 同：./run_long_lmp.sh ...
# 依次执行：捕获 → 原始/压缩/trimonly CCDG → SimGrid 验证 → 需求提炼+排布 →
#         无同步 CCDG 生成 → BookSim 两档注入（free 基线 + demand）
# 产物：runs/pipeline/short_lialocl_2688a_16r_<ts>/
#   ├── trace_16ranks_global.ccdg       (原始，含 setup)
#   ├── compact_16ranks_global.ccdg     (裁剪+折叠，供 SimGrid 验证)
#   ├── trimonly_16ranks_global.ccdg    (裁剪不折叠，需求编译器输入)
#   ├── validation_result.txt           (SimGrid 单步验证)
#   ├── demand_plan.log                 (需求提炼+排布报告)
#   ├── booksim_free/                   (free 基线档：trimonly 1 步)
#   └── booksim_demand/                 (demand 档：_demand.ccdg + _demand.est)
```

SimGrid 验证（validation_result.txt）在此工作流中的角色：**数据质量闸门**——它校验 COMPUTE 时长建模与消息集合守恒（迭代校准频率 ≤4 GHz 物理上限、T_simgrid 对 Loop time 误差 ≤5%），不是同步时序。demand 编译器主动丢弃时序/同步，但消费 compute_cycles/comm_bytes，提取失真会被 planner 静默吸收导致时间尺度整体失真，故进入排布前必须 PASS。

256r short 复用现有产物目录 `runs/pipeline/short_lialocl_2688a_256r_20260827_153703`（重跑或重新捕获亦可）。

## 4. 阶段 2：trimonly CCDG 生成

已内嵌为流水线 ③c（无条件生成，sim=任意值都产出）。手工重跑单阶段时参考（注意路径在运行目录下、不设 CCDG_COMPACT）：

```bash
cd /work1/jiangtao/lammps_trace
RUN_DIR=runs/pipeline/short_lialocl_2688a_16r_<ts>
env -u CCDG_COMPACT CCDG_TRIM_SETUP=1 \
  LD_LIBRARY_PATH=/work1/jiangtao/lammps_trace/install/lib \
  ./dumpi2ccdg/dumpi2ccdg "$RUN_DIR" \
    > "$RUN_DIR/trimonly_16ranks_global.ccdg" 2> "$RUN_DIR/ccdg_trimonly_gen.log"
```

**关键坑**：`dumpi2ccdg` 的环境变量用 `getenv()` **只查存在性不查值**：
- 要 trimonly 必须**不设置** `CCDG_COMPACT`（`CCDG_COMPACT=0` 仍会触发突发折叠），用 `env -u CCDG_COMPACT` 清除外层残留；
- 需要哪个模式就只设置哪个变量，切换模式用 `unset` 而非置 0。

## 5. 阶段 3：需求提炼（extract_demand）

已内嵌为流水线 ⑥（sim=demand 时）。手工重跑单阶段时参考：

```bash
cd /work1/jiangtao/lammps_trace/booksim2
python3 ccdg_demand.py runs/pipeline/short_lialocl_2688a_16r_<ts>/trimonly_16ranks_global.ccdg \
    -o runs/pipeline/short_lialocl_2688a_16r_<ts>/booksim_demand/demand_16ranks \
    --cap 2.5e10 --noc-ghz 2.0 --cpu-ghz 2.0 --hop-stride 4.0
```

提炼规则（`SYNC_TYPES`）：
- **丢弃**：WAIT / WAITALL / WAITANY / RECV / IRECV / BARRIER；
- **保留**：COMPUTE（compute_ops/compute_cycles）、SEND（comm_src/comm_dst/comm_bytes）、collective（ALLREDUCE/BCAST 等）。

守恒判据（验收标准 1）：
- `Σcomm_bytes(提炼后) == Σcomm_bytes(trace)` 精确相等；
- 方向直方图 ≤8 类（256r 三对主导：`(±1,0)/(0,±4)/(±8,0)`——LAMMPS 3D 分解在 16×16 mesh 上的投影）。

## 6. 阶段 4：时槽排布（plan_flow）

每 rank 两条虚拟时间线：

| 时间线 | 含义 |
|---|---|
| `pe_t[r]` | PE 链进度：COMPUTE 驻留、消息节点 +1 cycle |
| `port_t[r]` | 注入端口串行化：前一包尾 flit 时刻（端口 1 flit/cycle） |

排布规则：
1. **COMPUTE**：按 `node_dur` 推进 `pe_t`；
2. **SEND**：释放时刻 `t = max(pe_t, port_t)`，在 `PipelineLinkTable.first_free` 找 DOR 路径全程空闲的最早注入时刻，`reserve` 占用区间；
3. **collective（ALLREDUCE 等）**：BookSim 到达即展开、端口**连续 back-to-back** 注入全部轮次包 → 按 burst 整体平移搜索最早连续起点（逐轮独立排会低估注入反压）；
4. 全 rank round-robin 单节点推进，保持公平。

**PipelineLinkTable**（`--hop-stride 4.0`）：建模 BookSim IQ router pipeline（routing_delay=1 + vc_alloc_delay=1 + sw_alloc_delay=1 + traverse=1）。第 j 跳链路占用区间为 `[t + j×stride, t + j×stride + flits − 1]`。不用 stride 建模会导致真实路由器排队（congestion 上升）。

### 6.1 b 参数化相位注入许可（--phase，默认关）

Cerebras 论文机制：路由器"only when it is in the head state"才接受本地核注入，b+1 相位把 worker grid 分成 (b+1)×ny strips、每相位每 strip 一行一个核多播——**b 显式决定注入许可节奏**。本模式将其 P2P 同构加入 plan_flow：

| 概念 | 论文（多播） | 本实现（P2P 单播） |
|---|---|---|
| b | 多播传播距离 | 方向 d 的 DOR hop 数 `b_d = |dx|+|dy|` |
| 相位宽度 | 足迹 b+1 tiles | `w_d = b_d + 1`（DOR 路径占据 w_d 个连续坐标） |
| strips | 垂直条带宽 b+1 | 相位坐标 `c` = **垂直于传播轴的坐标**（水平方向用 x、垂直方向用 y）mod w_d |
| 注入许可 | HEAD 状态轮转 | `t ≡ c (mod w_d)` 格点；SEND 的 t0 = 最近许可格点 ≥ max(pe_t, port_t)；collective burst 起点对齐主导方向 w_d |

**零争用定理（同相位空间分离）**：同方向、同相位（c 相同）的两条消息，DOR 路径在空间上相距 ≥ w_d ⇒ 链路集合不相交 ⇒ 零冲突**与 flit 数/时间无关**。轴对齐方向严格成立；对角方向（corner 交换）同理由 x 分离保证。

**断言实现**：`PipelineLinkTable` 的每个已占用区间带 owner `(方向, 相位)` 标签，冲突时归因——`same_phase_conflicts` 必须为 0（理论断言），`cross_phase_conflicts`（相位间/异方向冲突）是交错格点的固有仲裁成本，由时槽表消化。

**与论文的语义差异**：论文的 b+1 相位是 **stage 内串行**（相位 p 完成后才轮到 p+1）；本实现的格点是**交错复用**（各相位每 w_d cycles 均可注入），相位间冲突靠时槽表仲裁。交错格点保留注入节奏确定性而不引入相位间 barrier，代价即 cross_phase_conflicts。

## 7. 阶段 5：无同步 CCDG 生成（emit）

产物：
- `<prefix>_demand.ccdg`：每 rank 链 `COMPUTE → SEND × N → collective × N`，`cross_rank_edges = []`，无任何同步节点，每节点带 `sched_est`；
- `<prefix>_demand.est`：`"node_id 释放时刻"` 两列表，被 BookSim `ccdg_schedule_file` 读取。

BookSim 零改动兼容性依据：
- SEND 无 `src_edge_indices` → 照常 `_injectPacket`（包到达时找不到 msg_id 仅打 WARNING，不阻塞）；
- WAIT/RECV 无 `dep_edge_indices` → 立即通过；
- ALLREDUCE 非阻塞展开（到达即注入、立即推进）；
- `sched_est` 门控每个节点的释放时刻（PE_GATED 分支）。

## 8. 阶段 6：BookSim 注入与结果

已内嵌为流水线 ⑥（sim=demand 时两档注入：free 基线 + demand）。手工重跑单档时参考：

```bash
cd /work1/jiangtao/lammps_trace/booksim2
RUN_DIR=/work1/jiangtao/lammps_trace/runs/pipeline/short_lialocl_2688a_16r_<ts>
CCDG_SCHED_FILE=$RUN_DIR/booksim_demand/demand_16ranks_demand.est \
./run_ccdg_mesh.sh $RUN_DIR/booksim_demand/demand_16ranks_demand.ccdg 3600 1
```

**命名规则**：`TAG = DIRTAG + _cap + _qd + _sched + _<steps>s + _<gating>`，其中 `DIRTAG = basename(dirname(ccdg))`。

**关键坑**：`DIRTAG` 取**父目录名**而非文件名——同一目录下放多个 CCDG 产物 TAG 相同、stats/CSV 互相覆盖。流水线已用 `booksim_free/`、`booksim_demand/` 两个子目录隔离两档；手工注入时每个 CCDG 用独立目录（运行目录名唯一，天然满足）。

结果读取：
- stats：`results/ccdg_mesh_<TAG>_stats.txt`（`total_sim_cycles` / `compute_cycles` / `blocked_cycles` / `congestion_cycles` / `sched_wait_cycles`）；
- CSV：`results/ccdg_mesh_results.csv`；
- 日志：`results/ccdg_mesh_<TAG>.log`。

账本守恒（验收标准 2）：

```
ranks × makespan = compute + blocked + congestion + sched_wait + PE_DONE_idle
```

- free 档：账本缺口 ≈ 0（16r 0.2%、256r 1.1%）；
- demand 档：`blocked = 0`，等待显性化为 `sched_wait`；残余 `congestion` 来自 collective 展开轮间注入反压（端口带宽成本，非同步）；账本缺口 = **PE_DONE idle**（rank 完成不均衡的自然闲置，无 BARRIER 会合的形态）。

## 9. 实测结果参考（2026-08-29，cap2.5e10，全部 unresolved=0、sent==recv）

| 档 | trace free | demand（本编译器） | 收缩 |
|---|---:|---:|---:|
| 16r short 4×4 | 225,887（blocked 1.96M，54.5%） | **201,878（blocked=0）** | **−10.6%** |
| 256r short 16×16 | 159,449（blocked 27.5M，68.1%） | **136,384（blocked=0）** | **−14.5%** |

planner 模型预测 256r = 136,382.4 vs 实测 136,384（误差 <0.01%）。

### 9.1 相位格点化三档验证（2026-08-29，同源载体，cap2.5e10，unresolved=0）

| 档 | trace free | demand 自由时槽 | demand 相位格点（--phase） | 格点化成本 |
|---|---:|---:|---:|---:|
| 16r short 4×4（trimonly 载体） | 229,926（blocked 56.1%） | 196,734（blocked=0） | **196,752（blocked=0）** | **+0.009%** |
| 256r short 16×16（compact 载体 153703，13824 包） | 180,669（blocked 52.1%） | 191,624（blocked=0） | **191,831（blocked=0）** | **+0.11%** |

- **零争用断言全 OK**：16r 10 方向、256r 14 方向（含对角 (8,3)/(8,±1)）`same_phase_conflicts` 全部 = 0——同相位空间分离定理在时槽表级严格成立；`cross_phase_conflicts`（相位间仲裁）16r ≈79、256r ≈3,229。
- 相位格点化成本 ≤0.11%（验收阈值 5%）——**确定性注入节奏几乎零代价**（自由时槽搜索本身已近似落在格点上）。
- est 预测 vs BookSim 实测偏差 ≤0.001%（188,131.8→196,752 的差值 4.6% 为网络实际路由开销，非排布误差）。
- **注意（载体差异）**：8/27 旧实验（§9 表，18368 包载体）demand 收缩 14.5%；本次同源补跑（13824 包新载体）256r demand 反超 free +6.1%——trace free 的天然错峰在部分载体上优于全局重排，需求驱动收益存在载体依赖边界，同源三档对比才是严谨口径。

| 档 | sched_wait | congestion | 账本缺口（PE_DONE idle） |
|---|---:|---:|---:|
| 16r demand | 1,166,691 | 351,408 | 15.8% |
| 256r demand | 17,165,144 | 370,439 | 25.0% |

## 10. CLI 参数参考

```bash
python3 ccdg_demand.py <trimonly.ccdg> \
    [-o 输出前缀（默认输入路径去扩展名）] \
    [--cap 2.5e10]      # 计算能力 ops/s（rate = cap / noc_ghz / 1e9）
    [--noc-ghz 2.0]     # NoC 频率 GHz
    [--cpu-ghz 2.0]     # CPU 频率 GHz（freq_ratio = cpu/noc）
    [--hop-lat 1.0]     # 单跳延迟（预留，模型以 stride 为准）
    [--ser-lat 1.0]     # 串行化延迟 cyc/flit
    [--slack 1.0]       # 飞行时间松弛因子（报告用，释放时刻由时槽驱动）
    [--hop-stride 4.0]  # 路由器 pipeline 每跳 cycles（routing+VA+SA+traverse）
    [--phase]           # Cerebras b 参数化注入许可：方向 d 只在 t ≡ c (mod w_d)
                        # 格点注入（w_d = b_d+1，c = 源 rank 垂直于传播轴的坐标）
                        # 默认关 = 自由时槽搜索；开时输出 phase lattice 表与
                        # 零争用断言（同相位冲突必须 = 0）
    [--phase-offset 0]  # 相位格点起点偏移（mod w_d，风险回退旋钮）
```

与 `run_ccdg_mesh.sh` 的模型常量一致：`flit_size_bytes=1`、`ser=1 cyc/flit`、频率 2.0 GHz（freq_ratio=1）。

## 11. v2 语义修正：blocked=0 从"删除语义"到"满足语义"（2026-08-31）

旧口径（§1–§10）删除 WAIT/RECV/BARRIER 后宣称 blocked=0 by construction——实测证明这是**删除语义**：跨 rank happens-before 随同步节点一起丢失，16r 载体 166 条依赖中 114 条（68%）存在"消费节点早于数据到达"（最严重差 5 万 cycle），free 档 52% 的 blocked 正是这些等待。v2 修复：

1. **happens-before 约束（默认开启）**：`extract_demand` 从 cross_rank_edges 提炼 `deps=[(send_id, consumer_id)]`（consumer = 接收 rank 上 WAIT 之后第一个保留节点）；`plan_flow` 对 consumer 施加 `est ≥ feeder_est + hops×stride + flits`，feeder 未规划时 defer（依赖图是真实执行 happens-before 的子图，无环）。plan 后全量复检断言 violations=0。`--no-hb` 保留旧乐观下界档。
2. **ejection 端口建模**：目的端口作为伪链路 `("ej", dst)`（offset = hops×stride）进时槽表，同节点同时到达互斥。仅 free-slot 模式；`--phase` 保持原断言域。
3. **cerebras wavelet 模式（`--cerebras`）**：`reorder_chain` 把每 rank 链按段（collective 界）重排为 X 阶段（dx≠0）先、Y 阶段后，段内按 (垂直坐标 c, dst) 排序；`plan_flow` 施加 X→Y 全局 stage barrier + b+1 相位串行轮转（tick q 在 tick q-1 全部落地后开放）。`--stage-barrier=global|rank|none`、`--wavelet-mode=serial|interleave`。
4. **pop-scan bypass + 链序一致**：barrier 下 Y send 会队头阻塞同 rank 后续 X send，循环改为"扫描第一个 ready 节点并 pop"；plan_flow 内部复制 chains，emit 前按 plan_t 重排链使 est 沿链单调（否则 BookSim 链序执行时 est 不可达，16r 曾差 83%）。
5. **载体 Z 压缩**：流水线捕获用 `processors K K 1`（Z 邻居本地化），方向集从 14 类（含 z 投影 (8,3) 类）缩到 8 类（面邻居 ±1 + PBC 接缝 ±(k-1)）。接缝是 torus→mesh 映射的真实代价（16r b=3、256r b=15），不是残留。

**Z 压缩载体四档实测**（lialocl 2688a，2026-08-31，全部 unresolved=0、sent==recv、断言全绿）：

| 档 | 16r | 256r | est↔实测 |
|---|---:|---:|---|
| free | 210,129 | 147,872 | — |
| hb | 193,572（+7.9%） | 137,010（+7.4%） | ≤0.02% |
| cerebras | 196,683（hb+1.6%） | 156,334（hb+14.1%） | ≤0.02% |
| ilv | 196,695（+0.01%） | 156,211（+0.08%） | ≤0.02% |

结论：需求驱动净收益（正确性口径）~+7.5% 且两规模一致；Cerebras 全机制在干净 2D 载体上 16r 近零代价、256r 代价主要来自接缝长波前的 stage 串行；相位串行 vs 交错几乎无差（代价在 barrier 不在轮转）。

## 12. fold-PBC（paper III-E 交叉排布）：接缝消解与多播依赖（2026-08-31）

论文 III-E 处理周期边界的方式：把周期维坐标环从两处剪开、两条半环**正序/倒序交错**铺上 fabric——环邻接全部变成 1–2 跳（接缝 1 跳），并断言 V-F "带宽不是限制资源，PBC 开关耗时相同"。实现为 `ccdg_demand.py --fold-pbc`（demand 档编译前 rank 置换 logical→fabric，节点 id/依赖不变，BookSim 零改动）与 `--fold-only IN OUT`（free 档 fold 基线重标号）。fold 后断言 max_b ≤ 2（16r: b1=64+32→混布，256r: b1=1248/b2=8784/max=2）。

**五档实测**（lialocl 2688a Z 压缩载体，全部 unresolved=0、断言全绿）：

| 档 | 16r | 256r |
|---|---:|---:|
| free | 175,046 | 151,285 |
| free_fold | 174,119（+0.5%） | 146,286（**+3.3%**） |
| hb（fold） | 163,844 | 137,867 |
| cerebras（fold） | 165,892（hb+1.25%） | 162,694（hb**+18.0%**） |
| ilv（fold） | 165,903 | 162,729 |

**反向发现：单播投影下 fold 对 cerebras 档净亏**（fold 前 +14.1% → fold 后 +18.0%）：

1. free 档 fold 赚——反应式调度只关心飞行时间，接缝 b=15（flight 64 cyc）消失是纯收益；
2. hb 档 fold 中性——接缝长飞行消失 ≈ 主流量 b=1→2 的占用增加，两者相抵；
3. cerebras 档 fold 亏——stage 串行把链路占用成本放大：fold 惠及的接缝消息只占 6%（624/10032），全体主流量（94%）b 翻倍则被串行调度乘上相位等待。

**与论文 "PBC is free" 的分歧定位**：论文的邻居变 2 跳发生在**多播波前**里（b=2 波前一次注入 fan-out 覆盖 3 tiles，不是 2 条单播），且 WSE 带宽富余；本复现是单播投影，b=2 = 实打实 2 链路 × stride 窗口。**推论：fold 的收益依赖多播原语支撑——MCAST 需求类型（§计划中）是复现论文 PBC 结论的前提，单播 fold 结论应标注为"多播缺失下的下界"。**

跨规模结论更新：Cerebras stage 串行代价 16r +1.3%、256r +14~18%，与 fold 与否弱相关——主因是串行化放大链路占用，不是接缝波前本身；相位串行 vs 交错在所有载体/模式上 ≤0.08%。
