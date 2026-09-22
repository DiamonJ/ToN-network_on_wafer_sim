# LiAlOCl 2688 原子强扩展计算—通信结果

## 1. 实验口径

- 场景：LiAlOCl 2688 原子，short-range，一次稳态 timestep。
- 拓扑：rank 与二维 mesh 节点一一映射。
- BookSim 模式：WSE-fast、fold-PBC on。
- 计算能力：25,000,000,000 ops/s。
- 静态计算量与通信量来自源码推导的无校准上下界。
- 硬件计算量采用 100 steps 的 perf `run(N)-run(0)`；本表为 1 repeat。
- 16/64 ranks 未超卖；256 ranks 在 240 个逻辑 CPU 上轻度 oversubscribe。

## 2. 核心结果

关键路径计算通信比定义为 `expected compute cycles / expected communication cycles`。

| 场景规模 | 关键路径计算/通信比 | 预期 cycles | BookSim cycles | 平均 packet 排队 | 平均 flit 排队 | 平均注入率 | 注入饱和占比 | 饱和时注入率 | BookSim 通信/计算比 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2688 atoms / 16 ranks (4x4) | 1.5083 | 124,467 | 124,474 | 0.000 cyc | 1293.000 cyc | 0.222 flit/slot | 22.203% | 1.000 flit/slot | 0.6259 |
| 2688 atoms / 64 ranks (8x8) | 0.4326 | 80,696 | 80,703 | 0.000 cyc | 1176.480 cyc | 0.305 flit/slot | 30.535% | 1.000 flit/slot | 2.7649 |
| 2688 atoms / 256 ranks (16x16) | 0.1142 | 84,897 | 84,904 | 0.000 cyc | 1141.800 cyc | 0.259 flit/slot | 25.853% | 1.000 flit/slot | 13.8519 |

## 3. 估算与实测计算量

理论 C1 是源码逐项推导的 algorithmic scalar-equivalent ops 区间；
实测 DP ops 来自 lane-weighted FP_ARITH_INST_RETIRED；退休指令包含
LAMMPS、MPI 和运行时执行的全部指令。三者不使用经验系数互相换算。

| Ranks | 理论总 C1 ops 区间 | 实测 DP ops/step | 实测总 instructions/step | 理论关键 rank ops 上界 | 预期 compute/comm cycles |
|---:|---:|---:|---:|---:|---:|
| 16 | 7,208,461–14,620,588 | 11,643,815 | 98,151,790 | 935,544 | 74,844/49,623 |
| 64 | 7,208,461–14,620,588 | 11,805,134 | 235,235,067 | 304,596 | 24,368/56,328 |
| 256 | 7,208,461–14,620,588 | 14,315,483 | 2,159,054,344 | 108,784 | 8,703/76,194 |

## 4. BookSim 指标定义

- **平均 packet 排队时间**：整包进入源端注入队列，到 head flit 成功注入网络的平均周期。
- **平均 flit 排队时间**：每个 flit 从进入注入队列到实际注入的平均周期，包含包内串行化等待。
- **平均注入率**：`injected_flits / (nodes × subnetworks × simulation_cycles)`。
- **注入饱和占比**：`backlogged_injection_slots / total_injection_slots`，表示注入端有积压的时间比例。
- **饱和时注入率**：`injected_flits / backlogged_injection_slots`，即有积压时每槽位实际发射率；它不是 synthetic offered-load sweep 的饱和拐点。
- **BookSim 通信/计算比**：全 rank 累计 `(blocked + congestion + sched_wait) / compute`。
- **关键路径计算/通信比**：编译器 `compute barrier / communication window`。它与上一项聚合方式不同，不能互为倒数。

## 5. 结果解读

1. 16 ranks 时关键路径计算/通信比为 1.5083，计算仍略占主导。
2. 64 ranks 时该比值降至 0.4326，通信开始主导；预期总 cycles 从 124,467 降至 80,696。
3. 256 ranks 时该比值仅 0.1142，通信窗口达到 76,194 cycles，预期总 cycles 回升到 84,897，说明 64–256 ranks 之间已经越过强扩展最佳点。
4. 三档实测 DP ops 为约 11.64M、11.81M、14.32M，与固定问题规模下理论总 C1 区间保持同量级。
5. 总退休指令从 98.15M 增至 2.159B，说明 rank 增加后 MPI/runtime 固定开销显著放大；256-rank 数据还包含轻度 oversubscription 的调度影响。
6. WSE 编排下 packet 排队为 0，饱和时注入率为 1 flit/slot；性能恶化主要来自通信窗口和全 rank 调度等待增长，而非源端注入队列无法服务。
   packet 排队为 0 表示 head flit 到队即发；flit 排队仍约 1,100–1,300 cycles，是长包尾部在单注入端口上的串行化等待。

## 6. 可复现性与限制

- 本次硬件 profile 使用 1 repeat，用于扩展趋势验证；正式统计建议设置 `COMPUTE_REPEATS=3` 后取逐 rank 中位数。
- 预期 cycles 使用理论 C1 上界作为保守计算预算，不使用硬件 profile 反推系数。
- 当前结果针对 short-range WSE 编排；PPPM/Kspace 长程通信需另建包含 FFT remap 的扩展性矩阵。
