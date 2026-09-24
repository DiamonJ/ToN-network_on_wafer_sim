# WSE Phase 2 v1 实现与验收记录

## 1. 阶段目标

本阶段完成 `实验平台重构方案.md` 中 Phase 2 的 v1：

> traffic-manager 级多播。Compiler 给出确定性的 wavefront 和共享链路足迹，
> manager 将足迹展开为 BookSim 可执行的 branch packet；BookSim 负责真实的
> flit、VC、router pipeline、credit 和 ejection 回放。

本阶段不实现 router 内部动态 flit 复制。后者属于可选的 Phase 2 v2。

## 2. 修改前的状态

Phase 1 完成时已经存在：

- `sim_type=wse`；
- `WSETrafficManager` 类；
- `.replay.ccdg/.replay.est`；
- compiler cycles 与 BookSim cycles 的 1% 误差检查。

但当时仍是兼容性桥接：

1. 每条原始 logical flow 被重新转换成端到端单播；
2. `WSETrafficManager` 只是 `CCDGTrafficManager` 的空子类；
3. manager 没有统计 wavefront、command 或 branch；
4. 没有验证共享 footprint edge 是否只回放一次；
5. 没有 congestion 验收；
6. 使用基础 cfg 的 8-byte flit 和 24 VC IQ profile。

因此它可以验证 EST，但还不等价于方案中定义的 manager-level footprint replay。

## 3. Compiler 后端修改

修改文件：

```text
wse_compiler.py
tests/test_wse_compiler.py
```

### 3.1 从 flow replay 改为 footprint replay

原实现：

```text
每条 logical flow
    → 一个端到端 BookSim SEND
```

新实现：

```text
每个 wavefront
    → 一个 command packet
    → 每个唯一 (directed link, depth)
        → 一个一跳 branch packet
```

Compiler 读取 `footprint_slots`，只对以下格式生成 branch：

```text
x0,y0 -> x1,y1
```

以下资源只参与编译期预留，不生成数据 branch：

```text
INJECT@x,y
EJECT@x,y
```

这样，一个被多个 logical flow 共享的 `(link,depth)` 在 BookSim 后端只生成一个
packet，而不是被多个端到端单播重复占用。

### 3.2 无 payload lineage 时的字节口径

当前 LAMMPS WSE Plan 只包含聚合消息大小，没有 payload/root lineage。Compiler
不能证明两个重叠 flow 携带完全相同的数据。

因此采用保守规则：

```text
同一个 (link,depth)：
    packet_count = 1
    packet_bytes = 所有经过该 slot 的 payload bytes 之和
```

这保证：

- 共享链路只创建一个 branch packet；
- 不会错误丢失独立 payload；
- logical message count/bytes 仍然守恒；
- link-level byte 数可以高于 logical payload byte 数，因为一个 payload 可能经过
  多个物理 hops。

256-rank case 中：

```text
logical messages       = 10,032
logical payload bytes  = 11,974,656
footprint branches     = 3,840
branch link bytes      = 22,477,056
```

### 3.3 Command wavelet

每个 wavefront 生成一个 command packet：

- payload bytes 为 0；
- BookSim 中仍占一个 flit；
- source 和 destination 都是 wavefront root node；
- 用于建模本地控制 wavelet，而不占用数据 mesh link；
- command 的注入和送达被单独统计。

### 3.4 Replay 资源闭包

仅有 tree-link 冲突闭包还不足以保证 BookSim 无 congestion，因为 manager 展开的
branch packet 仍需经过 PE injection 和 destination ejection。

因此 replay lowering 额外预留：

```text
INJECT@source
LINK@source->destination
EJECT@destination
```

对每个 packet 求不冲突的最早 release cycle。最终报告同时保留：

- `wse_semantic_total`：WSE wavefront 语义调度完成时间；
- `manager_replay_total`：加入 manager injection/ejection 约束后的时间；
- `compiled_total`：下游 BookSim 验收使用的完成时间。

### 3.5 WSE-fast profile

Compiler 新增：

```bash
--wse-fast-profile
```

它在用户指定的 `booksim2/ccdg_lammps_4x4.cfg` 基础上覆盖：

```text
flit_size_bytes = 4
routing_delay = 0
vc_alloc_delay = 1
sw_alloc_delay = 1
hop_stride_cycles = 3
```

BookSim IQ router 不允许 VC allocator 或 switch allocator 为 0 cycle。因此这里
使用当前 IQ 实现可执行的最浅 profile，没有伪造不可运行的单周期配置。

## 4. WSETrafficManager 修改

修改文件：

```text
booksim2/src/wse_trafficmanager.hpp
booksim2/src/wse_trafficmanager.cpp
booksim2/src/ccdg_trafficmanager.hpp
booksim2/src/ccdg_trafficmanager.cpp
booksim2/src/booksim_config.cpp
```

### 4.1 Replay 节点元数据

CCDG replay node 增加：

```text
wse_wavefront_idx
wse_kind = command | branch
wse_mode = multicast | reduction
```

普通 CCDG 节点使用默认值，因此原有 CCDGTrafficManager 行为不变。

### 4.2 生命周期 hooks

在 CCDG 网络执行引擎中增加两个虚函数：

```cpp
_OnPacketIssued(msg_id, node)
_OnPacketRetired(msg_id, destination)
```

普通 CCDG 实现为空操作；`WSETrafficManager` 覆盖它们进行 WSE 统计。

这样可以复用已经验证过的：

- packet/flit 创建；
- DOR routing；
- VC 分配；
- switch allocation；
- buffer 和 credit；
- tail retirement；
- network drain；

同时避免复制一份完整的 BookSim 网络循环。

### 4.3 Runtime 状态

`WSETrafficManager` 跟踪：

```text
expected_wavefronts
injected_wavefronts
msg_id → {wavefront, kind, mode}
commands_expected/delivered
branches_expected/delivered
multicast_branches_delivered
reduction_branches_delivered
branch_bytes_expected
```

packet 注入时：

1. 保存 `msg_id` 对应的 WSE 元数据；
2. 将 wavefront 加入 injected set。

tail flit 送达时：

1. 根据 `msg_id` 找到 WSE packet；
2. 增加 command 或 branch delivered；
3. 按 `multicast/reduction` 分类计数；
4. 删除已完成 packet 的运行时状态。

### 4.4 Congestion 指标

定义：

```text
wse_congestion_ratio =
    所有 PE_BACKPRESSURE cycles
    / (num_ranks × total_sim_cycles)
```

该指标验证 compiler 的 manager replay 闭包是否足以避免运行时注入拥塞。

### 4.5 新增配置入口

BookSim 配置增加：

```text
wse_program_file
```

当前执行仍使用 compiler 生成的 footprint replay carrier，但同时保留原始
`.program.json` 路径，用于标识 source WSE IR 和后续直接解析版本。

## 5. WSE runner 修改

修改文件：

```text
booksim2/run_wse_program.sh
```

Runner 现在会：

1. 读取 program 的 mesh、flit、compute capability 和 profile；
2. 设置 `sim_type=wse`；
3. 设置 `wse_program_file`；
4. 加载 footprint replay CCDG 和 replay EST；
5. 对 WSE-fast 设置 4 VC 和最浅 IQ pipeline；
6. 运行 BookSim；
7. 解析 WSE runtime stats；
8. 写入 `<prefix>.acceptance.json`；
9. 任一验收条件失败时返回非零状态。

验收条件：

```text
abs(booksim_cycles - compiler_cycles) / compiler_cycles <= 1%
wavefronts_injected == wavefronts_expected
branches_delivered == branches_expected
commands_delivered == commands_expected
wse_congestion_ratio <= 0.1%
```

## 6. Pipeline 修改

修改文件：

```text
run_noc_pipeline.sh
```

short pipeline 现在默认：

1. 捕获并合并 WSE Plan；
2. 生成静态 C1；
3. 使用 `--wse-fast-profile` 编译；
4. 运行 manager-level footprint replay；
5. 执行 cycles、wavefront、command、branch 和 congestion 闸门；
6. 将结果写入 `quality_gate.txt`。

失败信息不再只描述 EST 误差，而是明确覆盖：

```text
cycles / branch / command / congestion
```

long/PPPM 仍明确跳过该 short-range v1。

## 7. 新增统计

BookSim stats 文件新增：

```text
wse_wavefronts_expected
wse_wavefronts_injected
wse_commands_expected
wse_commands_delivered
wse_branches_expected
wse_branches_delivered
wse_multicast_branches_delivered
wse_reduction_branches_delivered
wse_branch_bytes_expected
wse_congestion_ratio
```

## 8. 16-rank 验收

Case：

```text
runs/pipeline/short_lialocl_2688a_16r_20260910_093912/
```

结果：

```text
compiler cycles        = 628,853
BookSim cycles         = 628,860
relative error         = 0.0011%
wavefronts             = 20 / 20
commands               = 20 / 20
branches               = 192 / 192
congestion ratio       = 0
status                 = PASS
```

## 9. 256-rank 验收

Case：

```text
cases/wse_short_256/
```

结果：

```text
compiler cycles                 = 143,547
BookSim cycles                  = 143,554
relative error                  = 0.0049%
wavefronts                      = 48 / 48
commands                        = 48 / 48
branches                        = 3,840 / 3,840
multicast branches delivered    = 1,920
reduction branches delivered    = 1,920
congestion ratio                = 0
status                          = PASS
```

验收详情：

```text
cases/wse_short_256/wse_phase1.acceptance.json
cases/wse_short_256/wse_phase1.booksim.stats
cases/wse_short_256/wse_phase1.report.json
```

## 10. 回归验证

完成以下检查：

- 3 个 Python compiler 单元测试通过；
- WSE-fast profile 参数测试通过；
- Python syntax check 通过；
- Bash syntax check 通过；
- BookSim C++ 增量编译成功；
- `make -q` 确认 build 最新；
- IDE linter 无新增错误；
- 16r 和 256r 端到端 v1 验收通过。

## 11. v1 的准确边界

本阶段已经实现：

- footprint 级共享 branch replay；
- command 单 flit建模；
- multicast/reduction branch 分类和守恒；
- BookSim 网络级 flit/VC/credit 回放；
- 无 congestion 的确定性 manager 调度。

本阶段没有实现：

- router 内部接收一个 flit并动态复制到两个 output；
- 分裂 credit 状态；
- router 内执行数值 reduction；
- payload/root lineage 驱动的真正相同 payload 去重；
- FFT、Grid3d 和 collective 的 WSE lowering。

因此当前结论是：

> Phase 2 v1 manager-level footprint replay 已完成并通过验收；Phase 2 v2
> router-level dynamic multicast/reduction 未实施。
