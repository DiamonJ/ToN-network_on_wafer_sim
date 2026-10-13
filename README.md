# ToN-network_on_wafer_sim

## 工作记录
### 20260910
1. lammps跑起来，了解运行lammps都需要哪些过程
2. 用之前发给你的in.lammps和data.lammps,跑一个case出来，了解in.lammps各行命令含义
3. 使用ai分析lammps计算流程和通信流程，对捕捉哪些计算过程和通信过程有一定了解，形成一个文档笔记
4. 学会使用dumpi或者其他profiling工具，对lammps运行过程中mpi流量进行捕捉，并设计一套方案对分子间作用力的计算量（flops/s）进行捕捉
5. 实现一个工具，可以直接从in.lammps配置文件，估算出每个rank的通信量T1和计算量C1
6. 通信量T1和计算量C1和profiling工具捕捉的T1'C1'进行比对验证，直到T1和C1收敛到合理误差范围内

**总结**：
已完成lammps的运行和数据采样，目前已有一版采样数据，并进行了校验（**校验过程待补充**）
### 20260917
1. 计算量改造: 通过捕捉指令数量，结合墙上时间和cpu频率估算计算量，同时对于通过代码理论计算方式给出计算通信量上下界
2. 在拓展性场景分析计算量变化，记录计算通信比和预期cycle数
3. 在booksim中记录平均排队时间，饱和注入率，通信计算时间比

### 20261009
1. 关注WSE gating的间隔完成相关排队时延实验
2. 测模拟原子的端到端时延
3. 尽可能扩大模拟规模，解决booksim大规模不稳定问题
4. 检查从原子到WSE通信策略链路是否集成到实验平台

#### BookSim 大规模稳定性与实验扩展

- **问题**：Linux 专用命令导致 macOS 假失败；高注入率不收敛或排空超时会中断整条扫描。
- **原因**：原运行器依赖 GNU `timeout` 和 `grep -P`，在 macOS 上会把成功仿真误报为失败；扫描脚本又把“不收敛、延迟超阈值、网络排空过久和真正执行错误”统一视为失败，并在单点失败后终止。大 mesh 每周期需遍历全部 router，排空成本随节点数和路径长度增长；若到达率超过网络服务率，队列还会持续增长，物理上不存在有限稳态时延。
- **处理**：用 Python 实现跨平台超时和进程组清理，以 `sed/awk` 解析结果；每个注入率点独立保存 JSON、支持断点续跑，并按 `PASS/SATURATED/TIMEOUT/NO_SAMPLES` 分类，单点失败不再影响后续扫描。绘图只保留第一个失败点之前的稳态结果并标出饱和率，同时分别记录源端排队时延和完整端到端 makespan，避免把超时或饱和瞬态当作稳态数据。
- **结果与边界**：真实 LAMMPS WSE plan 已验证到 256 ranks；基于实测通信模式的弱扩展投影已稳定完成 1024/4096 ranks，最大为 64×64 mesh、86016 原子、64512/64512 branches 送达。1024/4096 不是原生 trace，且尚未并行化 BookSim 核心推进算法。

#### 原子到 WSE 通信策略的集成情况

- 主入口 `run_noc_pipeline.sh` 默认开启 `WSE_PLAN_CAPTURE=1` 和 `WSE_COMPILE=1`：先由 `in.lammps` 生成原子体系并估算每-rank T1/C1，再在真实 LAMMPS 运行中通过 `LAMMPS_WSE_PLAN` 导出源码级通信计划，同时捕获 DUMPI trace。
- 每-rank WSE plan 会合并为 `wse_plan.json`，并与同源 trim-only CCDG 做逐方向消息数/字节数对拍；通过后由 `wse_compiler.py` 生成波前、EST 和 replay CCDG，再由 `run_wse_program.sh` 驱动 BookSim。
- BookSim 验收检查编译预测与实测 cycles 偏差 ≤1%、wavefront/branch/command 数量守恒和拥塞比例。现有 2688 原子、16-rank 样例中，WSE plan 与 CCDG 为 192 条消息、1175040 bytes 完全一致，编译/重放为 628853/628860 cycles，192/192 branches 送达。
- 当前 short-range 闭环已经默认集成并可自动失败退出；long/PPPM/Kspace 尚未进入 WSE 编译器。WSE 验收结果目前写入独立 acceptance 和 `quality_gate.txt`，但还没有作为独立档位并入最终统一 `evaluation.json`。

#### 实验图复现

```bash
python3 -m pip install matplotlib
python3 booksim2/sweep_synthetic_scaling.py --plot-only
```
