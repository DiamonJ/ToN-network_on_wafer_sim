# lmp-booksim 历史实验结果存档

本文件是 lmp-booksim skill 的历史结果数据存档。正文 SKILL.md 只保留最新口径的操作指引与结论；需要查具体数字、复核旧结论或对比新实验时才读本文件。

## 基线结果（2026-08，lialocl 2688a × 100 步，全部 unresolved=0）

| ranks | mesh | short cycles | long cycles | long/short |
|---|---|---:|---:|---:|
| 4 | 2x2 | 776,033,549 | 1,368,589,648 | 1.76× |
| 9 | 3x3 | 539,294,309 | 920,004,494 | 1.71× |
| 16 | 4x4 | 463,362,741 | 745,688,319 | 1.61× |
| 25 | 5x5 | 393,715,964 | 783,163,427 | **1.99×** |

结论：short cycle 单调下降（计算关键路径缩短）；long 在 16r→25r 反弹 +5.0%，Kspace 全局通信劣化吃掉扩展红利；long packets 16.4× 超线性增长（short 仅 8.3×）。

## 新流水线注入验证（2026-08-24，修复后 compact CCDG 直接注入，全部 unresolved=0）

| CCDG | mesh | steps | cycles/iter | 真实 Loop 折算(2GHz) | 一致性 |
|---|---|---:|---:|---:|---:|
| 16r short | 4x4 | 1 | 2,002,449 | ~0.96 ms | ✓ |
| 16r short | 4x4 | 10 | 1,981,913 | — | 单步偏差 -1.0% |
| 64r short | 8x8 | 1 | 1,061,308 | 0.53 ms | ✓✓ 几乎完全吻合 |
| 16r long | 4x4 | 1 | 7,744,011 | — | Kspace 注入正常 |

- **64r 8x8 mesh 的 BookSim cycles 与真实 Loop time 折算值几乎完全吻合**（0.531 vs 0.530 ms @2GHz）——强证据表明修复后 compact CCDG 在 mesh 上的计算+通信建模自洽，注入正常
- 新流水线 16r short cycles/iter ≈ 2.0M，旧 100 步基线 4.63M（高 2.3×）——旧基线受“setup 尾部间隙首 COMPUTE 未清零”失真影响，新值（≈真实时间）更可信
- long/short 比值 3.87×（16r 单步），方向与旧基线一致：Kspace 全局通信显著劣化

## BARRIER 锚定窗口修复验证（2026-08-25，compact_v2 直接注入，free gating 1 步，全部 unresolved=0）

| CCDG | mesh | 旧 cycles/iter | v2 cycles/iter | 变化 | P2P/rank | cycles/max_compute |
|---|---|---:|---:|---:|---:|---:|
| 4r short | 2x2 | — | 4,667,813 | — | 4（2×2，2 方向） | 1.042 |
| 16r short | 4x4 | 1,981,913 | 2,131,037 | +7.5% | 8（4×4，4 方向） | 1.088 |
| 64r short | 8x8 | 1,061,308 | 1,240,168 | **+16.8%** | 12（4×4×4 分解，6 方向） | 1.117 |
| 16r long | 4x4 | 7,724,138 | 8,476,967 | +9.7% | 38（含 kspace remap） | 1.205 |

- 64r +16.8% 为修复的直接证据：旧逻辑 halo P2P 全裁（0 SEND/rank），v2 恢复为每 rank 12 个 P2P 节点（6 方向 forward+reverse），NoC 争用被正确建模；此前“64r 完美线性/几乎完全吻合”结论因 P2P 被裁而部分失效，上方旧表数据均为旧逻辑产物
- cycles/max_compute 随规模单调上升（1.042→1.088→1.117→1.205），符合 mesh 跳数随节点数增长；16r long 最高（kspace FFT remap 通信密集）
- 计算量一致性：Σcompute/rank ÷ 2.5GHz ≈ Loop time（100–104%），SimGrid 频率校准仍落回物理范围
- 产物：`runs/pipeline/*/compact_v2_*.ccdg`（BARRIER 锚定版）；旧 `compact_*ranks_global.ccdg` 为 Loop-time 墙钟窗口版

## SimGrid 单步验证（2026-08-25，compact_v2 自动校准，全部 PASS）

| CCDG | T_real (Loop time) | 校准频率 | T_simgrid | E_total |
|---|---:|---:|---:|---:|
| 4r short | 1.782 ms | 2.515 GHz | 1.799 ms | 0.98% |
| 16r short | 0.771 ms | 2.595 GHz | 0.772 ms | 0.09% |
| 64r short | 0.424 ms | 2.701 GHz | 0.424 ms | 0.024% |
| 16r long | 2.800 ms | 2.512 GHz | 2.824 ms | 0.87% |

- **推翻旧结论“小步数验证必 FAIL”**：旧裁剪规则下单步误差 132×（setup 未裁净），BARRIER 锚定后单步验证 4/4 全部 ≤1%，校准频率落回物理范围（2.51–2.70 GHz）——单步 compact_v2 可直接作为 SimGrid 精度验证口径
- 用法：`python3 ccdg_simgrid_validate.py <run_dir> <num_ranks> <compact_v2_*.ccdg>`（自动校准含 SIMGRID_FREQ_CAP_GHZ=4.0 物理上限检查；固定频率用 SIMGRID_CPU_FREQ）

## BookSim free vs WSE gating（2026-08-25，compact_v2，1 步，全部 unresolved=0）

| CCDG | mesh | free | wse_p16_w2 | wse_p64_w2 | wse_p256_w2 | 最大增幅 | blocked@p256 |
|---|---|---:|---:|---:|---:|---:|---:|
| 4r short | 2x2 | 4,667,813 | 4,717,060 | 4,717,076 | 4,717,793 | +1.1% | 0.47M |
| 16r short | 4x4 | 2,131,037 | 2,178,488 | 2,189,544 | 2,204,510 | +3.4% | 1.77M |
| 64r short | 8x8 | 1,240,168 | 1,267,690 | 1,299,554 | 1,326,938 | **+7.0%** | 5.74M |
| 16r long | 4x4 | 8,476,967 | 8,622,951 | 8,806,698 | 8,861,034 | +4.5% | 12.3M |

- WSE 相位门控在本负载是**性能税**：cycles 增幅随 phase width 与规模增大（4r +1.1% → 64r p256 +7.0%），流量被挡在注入端拉长关键路径；blocked（累计）随规模超线性增长
- **旧 compact（P2P 裁掉）下 WSE 几乎无效果（64r blocked 仅 18–42K）——门控“无处发力”**；恢复 halo P2P 后 WSE 才真正起作用，旧“WSE 对 CCDG 负载无明显影响”结论需修正
- blocked 虽大（累计口径）但 cycles 增幅 ≤7%：阻塞多在非关键路径，关键路径受影响有限

## 需求驱动三档验证（2026-08-29，同源载体 1 步，cap2.5e10，全部 unresolved=0、sent==recv）

| 档 | trace free | demand 自由时槽 | demand 相位格点（--phase） | 格点化成本 |
|---|---:|---:|---:|---:|
| 16r short 4×4（trimonly） | 229,926（blocked 56.1%） | 196,734（blocked=0） | **196,752（blocked=0）** | **+0.009%** |
| 256r short 16×16（compact 153703，13824 包） | 180,669（blocked 52.1%） | 191,624（blocked=0） | **191,831（blocked=0）** | **+0.11%** |

- **零争用断言全 OK**：16r 10 方向、256r 14 方向（含对角 (8,3)/(8,±1)）same_phase_conflicts 全 0——同相位 DOR 路径空间相距 ≥w_d ⇒ 链路不相交，零冲突与 flit 数无关；相位格点化成本 ≤0.11%（验收阈值 5%），确定性注入节奏几乎零代价
- **收益存在载体依赖边界**：16r demand 收缩 −14.4%；但 256r 同源补跑中 demand 反超 free +6.1%（trace 天然错峰的关键路径短于全局重排的串行化空隙，contention slack = 链长的 129%）——旧 8/27 数据（18368 包另一载体）256r 收缩 −14.5%，同一规模不同捕获结论方向相反，**跨载体对比无意义**
- est 预测 vs BookSim 实测偏差 ≤0.001%；相位档 16r 196,752 / 256r 191,831 已入 CSV（booksim_phase 行）
