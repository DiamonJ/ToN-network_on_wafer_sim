# WSE Phase 1 实现与验收记录

## 1. 目标

本次工作按照 `实验平台重构方案.md` 完成 Phase 1 short-range WSE Compiler，
并补齐其验收所需的 BookSim manager-level 回放链路。

编译器不再使用 MPI trace 作为编译输入，而是读取：

1. LAMMPS 源码发射的 `wse_plan.json`；
2. `estimate_lammps_cost.py` 生成的 `static_cost_estimate.json`；
3. `booksim2/ccdg_lammps_4x4.cfg`。

最终目标是：

- 16-rank 和 256-rank short case 均通过全部编译断言；
- 编译器 EST 与下游 BookSim 实测完成时间的相对误差不超过 1%。

## 2. 实现内容

### 2.1 WSE Compiler

新增 `wse_compiler.py`，完成以下流程：

#### 输入规范化

- 解析 BookSim `key = value;` 配置格式；
- 校验 `topology=mesh`、`n=2`；
- 从 WSE Plan 的 `procgrid` 获取实际 worker grid；
- 将 cfg 中的 `k` 视为模板值，运行时 mesh 以 WSE Plan 为准；
- 从 cfg 读取 flit 大小、NoC 频率、计算能力和 router pipeline delay。

#### Placement 和 fold-PBC

- 要求处理器网格为 `Px × Py × 1`，完成 Z 压缩；
- 对偶数边长网格使用论文 III-E 的交叉 fold placement；
- 建立 LAMMPS logical rank 到物理 BookSim mesh node 的双射；
- 断言周期邻居在 fold 后不超过两个物理 hops。

#### H/V stage 编译

- `dimension=0` 编译为 H stage；
- `dimension=1` 编译为 V stage；
- 从 CommBrick round 几何得到 `b=max(round)+1`；
- 每个 stage 生成 `b+1` 个 strip phase；
- 按 perpendicular coordinate、phase tick 和方向组织 wavefront。

#### Multicast 和 reduction

- `forward`、`pair_forward` 编译为 multicast wavefront；
- `reverse`、`pair_reverse` 编译为 reduction wavefront；
- reverse wavefront 带有 `reduction_op=sum` 和分支节点信息；
- 每条 logical flow 保存 source、destination、bytes 和物理 DOR path。

#### TreeLinkTable

实现树状链路足迹的离线预留：

- 对 directed link、source injection 和 destination ejection 建立时隙；
- 使用 cfg 中的 router pipeline delay 计算 hop stride；
- 同一 `(link, depth)` 只生成一个 footprint slot；
- 因当前 WSE Plan 没有 payload lineage，重叠 slot 的负载取各 payload
  字节之和，避免错误地丢弃通信量；
- phase tick、H/V stage 之间使用 completion barrier 串联；
- 编译结束后检查所有方向链路无时间区间冲突。

#### Compute 模型

- 从静态估算结果读取每个 rank 的 `C1_steady_ops_midpoint`；
- 使用 cfg 的计算能力换算 NoC cycles；
- Phase 1 使用全局 compute barrier，再执行编译后的通信 stages；
- 细粒度 compute/communication overlap 尚未建模。

### 2.2 编译产物

对输出前缀 `<prefix>`，编译器生成：

- `<prefix>.program.json`：完整 WSE placement、compute、stage、wavefront、
  flow 和 footprint；
- `<prefix>.est`：`wavefront_id earliest_start_cycle`；
- `<prefix>.report.json`：统计、限制和编译断言；
- `<prefix>.replay.ccdg`：供 BookSim manager-level v1 使用的无同步 lowering；
- `<prefix>.replay.est`：replay node 的确定性 release time。

其中 replay CCDG 只是 WSE program 的后端载体，不是 MPI trace，也不重新引入
MPI WAIT、RECV 或 BARRIER 依赖。

### 2.3 BookSim 回放

新增：

- `booksim2/src/wse_trafficmanager.hpp`；
- `booksim2/src/wse_trafficmanager.cpp`；
- `booksim2/run_wse_program.sh`。

同时在 `TrafficManager::New` 中注册 `sim_type=wse`。

当前采用方案中推荐的 manager-level v1：

1. Compiler 将每个共享 `(link,depth)` footprint slot lower 为一次一跳
   branch packet，而不是把每个原始 flow 再做一次端到端单播；
2. 每个 wavefront 额外生成一个单 flit command packet；
3. `WSETrafficManager` 复用已经验证过的网络执行引擎，并独立统计
   wavefront、command、multicast/reduction branch；
4. BookSim 真实执行 DOR routing、IQ router pipeline、VC、credit、flit
   serialization 和 ejection；
5. replay 完成后，将实测 cycles 与 compiler report 中的
   `compiled_total` 比较；
6. cycles、计数守恒或 congestion 任一闸门失败时，runner 返回失败。

WSE-fast profile 使用 4-byte flit、4 VC，以及当前 BookSim IQ router 支持的
最浅流水：`routing=0`、`vc_alloc=1`、`sw_alloc=1`，对应 compiler
`hop_stride=3`。IQ router 不允许 VC/switch allocator 为 0，因此没有伪造
不可执行的单周期配置。

误差计算公式：

```text
relative_error = abs(booksim_cycles - compiler_cycles) / compiler_cycles
```

这不是路由器内动态 flit 复制版本。路由器级 multicast/reduction 属于可选的
Phase 2 v2。

### 2.4 Pipeline 集成

更新 `run_noc_pipeline.sh`：

- 新增 `WSE_COMPILE`，short 模式默认开启；
- WSE Plan 合并后自动运行 compiler；
- 自动调用 `run_wse_program.sh`；
- BookSim 误差超过 1% 时立即终止；
- `quality_gate.txt` 增加 `wse_compiler_booksim = PASS`；
- long/PPPM 明确跳过 Phase 1，而不是产生不受支持的结果。

## 3. 编译断言

`wse_phase1.report.json` 检查以下条件：

1. `placement_bijective`；
2. `z_compressed`；
3. `individual_routes_acyclic`；
4. `link_schedule_conflict_free`；
5. `est_monotonic`；
6. `message_count_conserved`；
7. `message_bytes_conserved`；
8. `fold_neighbor_hops_le_2`。

任一断言失败都会使编译命令失败。

## 4. 16-rank 验收

输入来自：

```text
runs/pipeline/short_lialocl_2688a_16r_20260910_093912/
```

编译结果：

- ranks：16；
- run messages：192；
- input/emitted bytes：1,175,040 / 1,175,040；
- stages：4；
- wavefronts：20；
- compiler cycles：628,853；
- BookSim cycles：628,860；
- relative error：0.0011%；
- wavefronts：20 / 20；
- branches：192 / 192；
- commands：20 / 20；
- congestion ratio：0；
- 所有八项编译断言：PASS。

验收文件：

```text
runs/pipeline/short_lialocl_2688a_16r_20260910_093912/
  wse_phase1.report.json
  wse_phase1.acceptance.json
```

## 5. 256-rank 验收

仓库中原有历史 256-rank BookSim 指标，但没有可复用的 256-rank
`wse_plan.json`。因此使用当前 LAMMPS WSE emitter 重新执行了真实的
`16 × 16 × 1` short LiAlOCl case，而不是复制或合成 16-rank 计划。

可复现输入：

```text
cases/wse_short_256/in.lammps
```

编译结果：

- ranks：256；
- run messages：10,032；
- input/emitted bytes：11,974,656 / 11,974,656；
- stages：4；
- wavefronts：48；
- compiler cycles：143,547；
- BookSim cycles：143,554；
- relative error：0.0049%；
- wavefronts：48 / 48；
- branches：3,840 / 3,840；
- commands：48 / 48；
- multicast/reduction branches：1,920 / 1,920；
- congestion ratio：0；
- 所有八项编译断言：PASS。

验收文件：

```text
cases/wse_short_256/
  wse_plan.json
  static_cost_estimate.json
  wse_phase1.program.json
  wse_phase1.report.json
  wse_phase1.acceptance.json
  wse_phase1.booksim.log
```

## 6. 验收结论

| Case | Compiler cycles | BookSim cycles | 相对误差 | 编译断言 |
|---|---:|---:|---:|---|
| 16r LiAlOCl short | 628,853 | 628,860 | 0.0011% | PASS |
| 256r LiAlOCl short | 143,547 | 143,554 | 0.0049% | PASS |

两组真实 short case 均满足：

```text
16r/256r 编译通过全部断言
BookSim relative error <= 1%
```

## 7. 使用方式

### 单独编译

```bash
python3 wse_compiler.py \
  <wse_plan.json> \
  <static_cost_estimate.json> \
  booksim2/ccdg_lammps_4x4.cfg \
  -o <output-prefix>
```

### BookSim 验收

```bash
booksim2/run_wse_program.sh \
  <output-prefix>.program.json \
  booksim2/ccdg_lammps_4x4.cfg
```

### 单元测试

```bash
python3 -m unittest tests.test_wse_compiler
```

### Pipeline

```bash
WSE_COMPILE=1 \
WSE_PLAN_CAPTURE=1 \
./run_noc_pipeline.sh short <ranks> <system> <atoms> <free|demand|both>
```

## 8. 验证工作

完成了以下检查：

- Python 单元测试通过；
- `wse_compiler.py` 和测试文件通过 Python syntax check；
- `run_noc_pipeline.sh`、`run_wse_program.sh` 通过 Bash syntax check；
- BookSim 增量构建成功且 build 状态最新；
- 修改文件无 IDE linter errors；
- 16r 和 256r BookSim 端到端验收均通过。

## 9. 当前边界

1. Phase 1 只支持 `forward/reverse/pair_forward/pair_reverse`；
2. FFT remap、Grid3d 和 collective 不在 short-range 论文机制范围内；
3. WSE Plan 暂无 payload/root lineage，因此相同链路上的独立 payload 仍需累加；
4. Compute 目前是静态全局 barrier，没有细粒度 overlap；
5. 当前 BookSim 后端是 manager-level footprint branch replay，不是
   router-level 动态 flit 复制；
6. `--wse-fast-profile` 在指定基础 cfg 上覆盖为 4-byte flit、4 VC runner
   配置和 IQ router 可支持的最浅流水。

这些限制均会在编译报告中显式记录，不影响本次 Phase 1 验收结论。
