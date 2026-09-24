***

# 🚀 TaskList: DUMPI + CCDG + BookSim 闭环仿真流水线

注意，所有新建的文件都必须放在/work1/jiangtao/lammps\_trace/目录下，不得干扰其他任务和路径

## 阶段一：环境搭建与 DUMPI Trace 提取 ✅

### Task 1.1: 编译与部署 DUMPI 基础库 ✅

- **目标**：从源码构建 `sst-dumpi` 动态链接库及 C 头文件。
- **输入**：`sst-dumpi` 源码仓库 (SST-Core / SST-Macro 官方组件)。
- **动作**：

1. 执行 `bootstrap.sh` 并配置 Makefile：`./configure --prefix=$INSTALL_DIR --enable-dumpi`
2. 编译并安装 `libdumpi.so` 和 `libdumpi.a`。
3. 确认生成头文件 `dumpi/libdumpi/dumpi.h` 与提取工具 `dumpi2ascii`。

- **产出**：动态库 `libdumpi.so` 及 C++ 头文件路径。
- **验证标准**：运行 `dumpi2ascii --help` 无报错；检查 `libdumpi.so` 符号表存在 `dumpi_start_stream_read`。

### Task 1.2: 运行 LAMMPS 并捕获 DUMPI Trace

- **目标**：无侵入式捕获 LAMMPS 运行轨迹。
- **输入**：标准 MPI 版 LAMMPS 可执行文件（如 `lmp_mpi`）、LAMMPS 输入脚本（如 `in.lj`，Rank 数设为 4 或 16）。
- **动作**：

1. 编写运行脚本，注入动态库环境变量，使用命令来调用lammps：root\@node77:/work1/jiangtao/DPMD-V1# sudo ./sw/script/run\_lammps\_test-jt.sh --lmp-case-dir ./vrf/sys\_vrf/PY\_REF\_V2/lammps\_config/ --lmp-mpi-ranks 16 --mdp-mode off --lmp-steps 100

   可能会涉及到run\_lammps\_test-jt.sh的修改

```bash
mpirun -np 16 LD_PRELOAD=/path/to/libdumpi.so ./lmp_mpi -in in.lj

```

1. 捕获每个 Rank 生成的 `.dumpi` 二进制文件与 `.meta` 元数据文件。

- **产出**：`dumpi-*.bin` 文件列表及 `.meta` 文件。
- **验证标准**：使用 `dumpi2ascii dumpi-*.meta` 能够成功打印出包含 `MPI_Send`、`MPI_Recv` 及 `dumpi_clock` 的文本日志。

***

## 阶段二：Trace 解析与依赖图（CCDG）转换工具开发

### Task 2.1: 编写基于 `libdumpi` 的 C++ 解析器框架

- **目标**：构建二进制 Trace 解析流，提取计算耗时与 MPI 通信事件。
- **输入**：`dumpi-*.bin` 及 `dumpi.h`。
- **动作**：

1. 创建 C++ 项目 `dumpi2ccdg`，链接 `libdumpi.a`。
2. 注册 DUMPI C 回调函数：`dumpi_on_send`、`dumpi_on_recv`、`dumpi_on_isend`、`dumpi_on_irecv`、`dumpi_on_wait`。
3. 提取 CPU 耗时计算逻辑：通过相邻两次 MPI 调用的 `cpu_time` 差值计算 `Compute_Cycles`。

- **产出**：`dumpi2ccdg` 可执行解析工具。
- **验证标准**：输入一个 `.meta` 文件，控制台能顺序输出每个 Rank 的事件序列（格式如：`[Rank 0] COMPUTE 12000 cycles -> ISEND to Rank 1 (4096 bytes, Tag 10)`）。

### Task 2.2: 实现 CCDG 依赖推导引擎

- **目标**：将事件流重构为有向无环图（DAG），关联 RAW（Read-After-Write）通信依赖。
- **输入**：Task 2.1 提取的事件流。
- **动作**：

1. 维护匹配数据结构：`PendingSendMap <(Src, Dst, Tag), SendEvent>` 与 `PendingRecvMap <(Src, Dst, Tag), RecvEvent>`。
2. 识别非阻塞通信依赖：将 `MPI_Isend`/`MPI_Irecv` 的句柄与后续的 `MPI_Wait` / `MPI_Waitall` 进行绑定。
3. 导出标准 Intermediate Representation (IR) 格式（JSON 或二进制 `.ccdg`），结构包含：

- `NodeID`, `RankID`, `NodeType` (`COMPUTE` / `COMM`)
- `ComputeCycles`
- `CommSrc`, `CommDst`, `MsgSizeBytes`, `Tag`
- `Predecessors` (前驱依赖 NodeID 列表)
- **产出**：`.ccdg` 结构化依赖图文件。
- **验证标准**：通过 Graphviz 输出简单 2-Rank 交互的 `.dot` 图，人工核对 `MPI_Wait` 节点的前驱是否正确指向对应的通信与计算节点。

***

## 阶段三：BookSim 2.0 驱动层扩展 (`CCDGTrafficManager`)

### Task 3.1: 继承扩展 BookSim `TrafficManager`

- **目标**：在 BookSim 中实现基于 CCDG 依赖图的状态机驱动器。
- **输入**：BookSim 2.0 源码、Task 2.2 生成的 `.ccdg` 解析代码。
- **动作**：

1. 在 BookSim 源码 `src/` 下新增 `ccdg_trafficmanager.hpp` 和 `.cpp`。
2. 定义 PE（Processing Element）状态结构体 `PERankState`（状态：`COMPUTE`, `INJECTING`, `BLOCKED`）。
3. 实现配置项解析：引入 `frequency_ratio`（CPU 主频/NoC 主频）、`flit_size_bytes` 以及 `rank_mapping_file`。

- **产出**：包含 `CCDGTrafficManager` 类声明与配置项注册的代码。
- **验证标准**：修改 BookSim `booksim_config`，设置 `traffic = ccdg` 能成功被 BookSim 识别并实例化。

### Task 3.2: 实现周期精确推进（`_Step`）与包注入逻辑

- **目标**：根据依赖图驱动 BookSim 主时钟周期推进。
- **输入**：`CCDGTrafficManager` 框架。
- **动作**：

1. 重写 `_Step()` 函数：

- 遍历所有 PE，若处于 `COMPUTE` 状态，对 `remaining_cycles` 递减。
- 当计算周期归零，读取 CCDG 中的下一个节点：
- 若为 `COMM` 注入节点，计算 $\text{FlitNum} = \lceil \text{MsgSizeBytes} / \text{FlitSize} \rceil$，调用 `_GeneratePacket()` 写入 BookSim 注入队列。
- 若为 `WAIT` 阻塞节点，检查依赖列表，若有未到达消息，设置 PE 为 `BLOCKED`。
- **产出**：带依赖控制逻辑的 `_Step()` 推进代码。
- **验证标准**：当 PE 处于 `BLOCKED` 时，即使过去多个 Clock Cycle，该 PE 也不会推进下一个计算节点。

### Task 3.3: 实现包到达回调（`_RetireFlit`）与依赖解锁

- **目标**：基于物理网络实际传输延迟解锁软件依赖。
- **输入**：BookSim 的 `_RetireFlit` 虚函数重写。
- **动作**：

1. 重写 `_RetireFlit(Flit *f, int dest)`：

- 检查 `f->tail` 是否为 `true`（消息尾包到达）。
- 从 Flit Header 中解析 `src`, `tag`, `msg_id`。
- 查找 `dest` 对应的 `PERankState`，将对应的依赖记录从 `pending_dependencies` 集合中清除。
- 若 `pending_dependencies` 为空，将该 PE 状态由 `BLOCKED` 恢复为 `COMPUTE`。
- **产出**：完整的闭环反压控制代码。
- **验证标准**：单步调试下，当 Tail Flit 穿过 NoC 到达 Ejection Port 的瞬间，目标 PE 状态立即更新为 `COMPUTE`。

***

## 阶段四：验证与端到端闭环测试

### Task 4.1: 单元验证（Ping-Pong & Ring Benchmark）

- **目标**：排除 LAMMPS 复杂逻辑干扰，验证驱动器逻辑正确性。
- **输入**：编译一个 2-Rank Ping-Pong MPI 小程序，提提取 `.dumpi` 并转为 `.ccdg`。
- **动作**：

1. 在 2x2 Mesh 拓扑的 BookSim 上运行该 CCDG。
2. 打印每个 Cycle 两个 Rank 的状态转换日记。

- **产出**：Ping-Pong 测试的时序 Log。
- **验证标准**：Rank 1 的接收开始时间严格等于 Rank 0 的 Flit 传输完成时间 + NoC 物理延时。

### Task 4.2: LAMMPS 真实负载端到端仿真

- **目标**：完成 LAMMPS 仿真并在 BookSim 中提取 NoC 物理性能指标。
- **输入**：Task 1.2 提取的 16-Rank LAMMPS CCDG 依赖图，配置 4x4 Mesh NoC 参数。
- **动作**：

1. 执行 BookSim 仿真：`./booksim config_lammps_4x4.cfg`。
2. 收集统计指标：总仿真 Cycle 数、平均 Flit 延迟、VC 占用率、Crossbar 冲突率。

- **产出**：全闭环仿真报告（包含硬件 NoC 瓶颈对 LAMMPS 运行时间的真实拖慢比例）。
- **验证标准**：仿真无死锁（Deadlock）顺畅运行至最后一个 Compute 节点，打印出完整的仿真周期统计。

