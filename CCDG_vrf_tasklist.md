# 注意，所有新建的文件都必须放在/work1/jiangtao/lammps\_trace/目录下，不得干扰其他任务和路径


SimGrid 的编译和安装建议优先使用**源码编译**（确保启用 `SMPI` 和 `smpi_replay` 功能），同时也提供了系统包管理器的快速安装备选路径。

---

## 包含安装步骤的完整 TaskList

```
[Task 0: SimGrid Setup] ──► [Task 1: Baseline Data] ──► [Task 2: DUMPI->CCDG Parser]
                                                                  │
[Task 5: Validation Gate] ◄── [Task 4: SimGrid Platform & Replay] ◄── [Task 3: CCDG->SimGrid Exporter]

```

1. **Task 0: SimGrid 环境准备与源码安装 (Environment Setup):** 目标：在本地/服务器部署包含 SMPI 模块的 SimGrid 工具链.
* **任务描述**：检查系统依赖，编译安装包含 SMPI 和 `smpi_replay` 的 SimGrid，并配置环境变量。
* **输入**：Linux 环境（Ubuntu 20.04/22.04/24.04 或 CentOS/RHEL）。
* **操作步骤**：
1. **安装构建依赖**：



```bash
sudo apt-get update && sudo apt-get install -y \
    git cmake build-essential libboost-dev libboost-tools-dev \
    libboost-thread-dev libboost-system-dev flex bison

```

2. **下载与编译 SimGrid 源码**：

```bash
git clone https://framagit.org/simgrid/simgrid.git simgrid-src
cd simgrid-src
mkdir build && cd build
cmake -DENABLE_SMPI=ON -DENABLE_DOCUMENTATION=OFF -DENABLE_MALLOCATOR=OFF ..
make -j$(nproc)
sudo make install
sudo ldconfig

```

*(备选方案：若系统环境受限，可尝试直接使用包管理器 `sudo apt-get install simgrid`)*

* **交付物 (Outputs)**：
* 系统 PATH 中可调用的 `smpicc`、`smpicxx`、`smpirun` 和 `smpi_replay` 二进制工具。


* **验收标准 (Acceptance Criteria)**：
* 执行 `smpirun --version` 输出正确版本号（通常为 3.32+ 或 Git 主干版本）。
* 执行 `smpi_replay --help` 返回 Exit Code 0，无动态库缺失报错 (`libsimgrid.so` 加载正常)。




2. **Task 1: 基准数据采集 (Baseline Data Acquisition):** 目标：获取 LAMMPS 真实运行耗时与 DUMPI 原始 Trace.
* **任务描述**：配置并运行真实的 LAMMPS 任务，同时挂载 DUMPI 库，记录真实物理耗时 $T_{\text{real}}$ 并导出二进制 Trace。
* **输入**：
* LAMMPS 模拟脚本 `in.lj`（如 1000 步 Lennard-Jones 体系）。
* DUMPI 拦截库 `libdumpi.so`。


* **操作步骤**：
1. 执行命令：`mpirun -np N LD_PRELOAD=libdumpi.so ./lmp_mpi -in in.lj`
2. 解析 LAMMPS 标准输出的 `Loop time` 得到 $T_{\text{real}}$。
3. 收集生成的 `dumpi-*.bin` 文件。


* **交付物 (Outputs)**：
* `dumpi_traces/` 目录（包含各 Rank 的二进制 trace）。
* `baseline_metrics.json`：包含 $T_{\text{real}}$、Rank 数量、总步数。


* **验收标准 (Acceptance Criteria)**：
* `baseline_metrics.json` 中的 $T_{\text{real}} > 0$ 且生成的 DUMPI 文件数量等于 Rank 数量 $N$。




3. **Task 2: 基于已有dumpi2ccdg模块，解析 DUMPI Trace，提取计算/通信依赖有向无环图.

* **输入**：`dumpi_traces/` 目录。
* **操作步骤**：
1. 解析点对点通信（`MPI_Isend`, `MPI_Irecv`, `MPI_Wait`）与集合通信（`MPI_Allreduce` 等）。
2. 计算邻近 MPI 调用之间的 CPU 耗时，建模为 `Compute` 节点。
3. 构建边依赖：单 Rank 序向依赖边 + 跨 Rank 消息匹配边（Matching `Send` to `Wait/Recv`）。
4. 拓扑排序校验，确保无环路。


* **交付物 (Outputs)**：
* `ccdg.json`（包含节点类型、关联 Rank、数据 Payload 大小、CPU 计算耗时、前驱/后继节点 ID）。


* **验收标准 (Acceptance Criteria)**：
* DAG 完整性检测：`ccdg.json` 能够通过 Kahn 拓扑排序算法（无环）。
* 点对点匹配率：100% 的 `Isend` 节点在图中有且仅有一个对应的 `Recv/Wait` 消费边。




4. **Task 3: CCDG 转 SimGrid smpi_replay 导出器:** 目标：将 CCDG 图导出为 SimGrid 原生 Replay 日志.
* **任务描述**：编写转换脚本，将 `ccdg.json` 导出为 SimGrid `smpi_replay` 工具支持的文本 Trace 格式（每 Rank 一个 action 文件）。
* **输入**：`ccdg.json`。
* **操作步骤**：
1. 映射 CCDG 节点到 SimGrid 动作语法：
* 计算节点 $\rightarrow$ `compute <seconds>`
* 发送节点 $\rightarrow$ `send <dst_rank> <size_bytes>`
* 接收节点 $\rightarrow$ `recv <src_rank>`
* 集体通信 $\rightarrow$ `allreduce <size_bytes>`


2. 按 Rank 拆分为 `actions_0.txt`, `actions_1.txt` ...
3. 生成包含文件映射关系的 `replay.meta`。


* **交付物 (Outputs)**：
* `simgrid_traces/` 目录（包含 `replay.meta` 及各 Rank 的 `actions_*.txt`）。


* **验收标准 (Acceptance Criteria)**：
* 导出的轨迹文件语法格式无误，总 `send` 动作次数与总 `recv` 动作次数严格相等。




5. **Task 4: SimGrid 平台建模与回放执行:** 目标：配置硬件环境并运行 SimGrid 回放得到 T_simgrid.
* **任务描述**：根据 Task 1 运行机器的物理参数（CPU 算力、网卡带宽、延迟）编写 SimGrid `platform.xml`，使用 `smpi_replay` 重放轨迹。
* **输入**：
* `simgrid_traces/` 目录。
* 宿主机硬件参数（如 100Gbps 带宽、500ns 延迟、算力 Scaling 系数）。


* **操作步骤**：
1. 生成 `platform.xml`（定义节点配置与 Link 参数）。
2. 执行回放命令：



```bash
smpirun -platform platform.xml -hostfile hostfile.txt \
        --cfg=smpi/replay-timing:yes smpi_replay replay.meta

```

3. 从 SimGrid 提取虚拟结束时间 $T_{\text{simgrid}}$。

* **交付物 (Outputs)**：
* `simgrid_result.json`（包含 $T_{\text{simgrid}}$、各节点计算耗时、通信耗时统计）。


* **验收标准 (Acceptance Criteria)**：
* `smpi_replay` 正常退出（Exit code 0），无死锁或未匹配通信的崩溃报错。




6. **Task 5: 自动化误差评估门控 (Validation Gate):** 目标：对比 T_real 与 T_simgrid，判断 CCDG 是否有效.
* **任务描述**：计算真实耗时与 SimGrid 仿真耗时的相对误差，判定 CCDG 依赖提取与计算时间量化的正确性。
* **输入**：`baseline_metrics.json` 和 `simgrid_result.json`。
* **操作步骤**：
1. 计算相对误差：

$$E_{\text{total}} = \frac{\vert{}T_{\text{real}} - T_{\text{simgrid}}\vert{}}{T_{\text{real}}} \times 100\%$$


2. 评估临界路径（Critical Path）上的计算/通信比例对齐度。
3. 输出诊断分析报告。


* **交付物 (Outputs)**：
* `validation_report.json`（包含 $E_{\text{total}}$ 及 Pass/Fail 结论）。


* **验收标准 (Acceptance Criteria)**：
* **PASS 门槛**：$E_{\text{total}} \le 5\%$（或工程允许的 $8\%$ 范围内），判定 CCDG 构建有效，解除阻断，允许进入阶段二（BookSim 中间件开发）。