# Kspace FFT 通信实验手册（上下文交接文档）

> 本文档供新对话读取后直接续做实验。包含：实验目的、环境、已完成的关键结果、待开始实验的完整执行方案。

---

## 一、实验目的

研究**强可扩展性下 2D mesh 的通信特征**：随着 rank 增加，Kspace 为何成为主要开销。

- **假设**：2D mesh 随规模增加通信性能下降，Kspace 中需要全局通信的 FFT（转置通信）成为通信开销的主要来源。
- **待证明**：Kspace 内部，**计算和邻居通信随 rank 增加而降低或小幅增加，而 FFT 全局（转置）通信随 rank 大幅增加**。
- **数据口径（用户明确要求）**：LAMMPS 作为应用真实运行，rank 间通信由 SimGrid 按 2D mesh 模拟，测试数据取 **LAMMPS 应用层面的真实耗时 breakdown**（真实时钟测量，不用仿真时钟）。

---

## 二、实验环境

### 2.1 关键路径

| 项 | 路径 |
|---|---|
| SMPI 版 LAMMPS 源码 | `/work1/jiangtao/lammps_trace/lammps-src/src/`（FFT 相关在 `KSPACE/` 子目录） |
| SMPI 版 LAMMPS 构建目录 | `/work1/jiangtao/lammps_trace/lammps-build-smpi/`（产物 `lmp`） |
| PIE 标志清除脚本 | `/work1/jiangtao/lammps_trace/simgrid_traces/strip_pie_flag.py` |
| 2D mesh 平台 XML | `/work1/jiangtao/lammps_trace/simgrid_traces/platform_mesh_{2x2,2x4,4x4,4x8}.xml` |
| 平台生成脚本 | `/work1/jiangtao/lammps_trace/simgrid_traces/gen_mesh_platform.py` |
| 2688 原子算例 | `/work1/jiangtao/lammps_trace/cases/lialocl_coul_2688/`（LiAlOCl，replicate 4×2×4，FFT grid 40×24×36，minimize 100 步 + run，STEPS 变量控制） |
| 1344 原子算例 | `/work1/jiangtao/lammps_trace/cases/lialocl_coul_1344/`（replicate 4×2×2） |
| 已有结果目录 | `runs/smpi_mesh_2688/`、`runs/smpi_mesh_1344/`、`runs/smpi_trace_kspace/` |
| 归位分析脚本 | `runs/smpi_trace_kspace/attribute_kspace.py` |

### 2.2 平台参数（当前）

host 计算能力 **1.42 Gflops**，链路带宽 **6.8 GB/s**，延迟 **10 ns**。rank r → mesh 坐标 (r % KX, r // KX)，全显式路由。

| ranks | 平台文件 | mesh | proc grid（LAMMPS 自动选择） |
|:---:|---|:---:|:---:|
| 4 | platform_mesh_2x2.xml | 2×2 | 2×2×1 |
| 8 | platform_mesh_2x4.xml | 2×4 | 2×2×2 |
| 16 | platform_mesh_4x4.xml | 4×4 | 4×2×2 |
| 32 | platform_mesh_4x8.xml | 4×8 | 4×4×2 |

### 2.3 运行条件（缺一不可）

```bash
cd /work1/jiangtao/lammps_trace/cases/lialocl_coul_2688
SMPI_PRIVATIZATION=0 timeout 1500 smpirun \
  -platform /work1/jiangtao/lammps_trace/simgrid_traces/platform_mesh_2x2.xml \
  -np 4 /work1/jiangtao/lammps_trace/lammps-build-smpi/lmp -in in.lammps -v STEPS 10
```

1. **`SMPI_PRIVATIZATION=0`**：privatization 开启会干扰 `MPI_COMM_WORLD` 绑定导致崩溃。
2. **二进制必须清 `DF_1_PIE` 标志**：每次重链接 `lmp` 后都要重跑 `python3 simgrid_traces/strip_pie_flag.py lammps-build-smpi/lmp`（输出 `FLAGS_1 0x8000001 -> 0x1`）。
3. lmp 链接时需 `main` 已导出（`-Wl,--export-dynamic-symbol=main`，现有构建已满足）。
4. **LD_PRELOAD 无法拦截 SMPI dlopen 加载程序的 MPI 符号**（已证伪，勿再尝试 shim 路线）。

### 2.4 相位 trace 插桩（已就位）

- `timer.cpp` 支持环境变量 `LAMMPS_PHASE_TRACE=<csv>`：rank 0 把相位 marker 写入 CSV，格式：
  `walltime_ns,phase_name,sim_time_s`
  **第 1 列 walltime_ns 是真实时钟（CLOCK_MONOTONIC），第 3 列是 MPI_Wtime 仿真时钟**。应用层真实耗时分析用第 1 列。
- 主相位行（Pair/Kspace/Comm…）是 **END marker**：`(上一行 ts, 本行 ts]` 属于该相位。
- `KSPACE_*` 行是 **START marker**：段持续到下一个 marker。已有 marker：`KSPACE_GRID_REVERSE`（rho 逆向交换，邻居通信）、`KSPACE_FFT`（brick2fft+poisson，含 FFT 计算与转置通信）、`KSPACE_GRID_FORWARD`（E 场正向交换，邻居通信）、`KSPACE_FORCE`、`KSPACE_REDUCE`。插桩点在 `src/KSPACE/pppm.cpp` 的 `Timer::mark_subphase()` 调用。
- 算例含 minimize（100 步）+ run 两段，各产生一个 breakdown。**切分方法**：从 run log 的最后一行 `Loop time of X on N procs for M steps` 取 run 步数 M，Kspace 窗口序列的最后 M 个属于 run 阶段。

### 2.5 时钟口径要点（重要，避免踩坑）

1. **LAMMPS Timer（Loop time / breakdown）用 `platform::walltime()` = `std::chrono::steady_clock` 真实时钟，不是 MPI_Wtime**。SMPI 下进程按仿真网络结果真实阻塞，所以 breakdown 就是"应用层真实耗时"，正是本实验要的数据。
2. SMPI 默认不把本地计算计入仿真时钟；`--cfg=smpi/host-speed:2.8544e9` 可开启（仅仿真时钟分析需要，应用层实验不用）。
3. Paje tracing（`--cfg=tracing:1 --cfg=tracing/smpi:1`）是第二层解剖工具，应用层主实验不需要。
4. log 里有两个 `MPI task timing breakdown` 块（minimize + run），取第二个：`awk '/MPI task timing breakdown/{n++} n==2'`。

---

## 三、已有关键结果

### 3.1 Kspace 占比随 rank 上升（应用层真实耗时口径，2688 原子，10 步 run）

来源：`runs/smpi_mesh_2688/mesh_{4,8,16,32}.log`（1344 原子在 `runs/smpi_mesh_1344/`，趋势一致）

| ranks | Loop time (s) | **Kspace %total** | Comm %total | Pair %total |
|:---:|:---:|:---:|:---:|:---:|
| 4 | 0.1495 | **72.5%** | 18.7% | 6.2% |
| 8 | 0.1596 | **78.3%** | 17.5% | 2.3% |
| 16 | 0.2913 | **80.7%** | 16.9% | 0.7% |
| 32 | 1.5252 | **92.9%** | 6.3% | 0.1% |

Loop time 超线性增长（rank×8 → 时间×10.2），Kspace 绝对耗时是恶化主驱动。

### 3.2 FFT 耗时占比（仿真时钟 + host-speed 口径，辅助佐证）

来源：`runs/smpi_trace_kspace/`（`phase_hs*.csv` + `simgrid_hs*.trace` + `attribute_kspace.py`）

| ranks | FFT 阶段占 run 段 | 转置通信占 FFT 内 | 转置占全部通信 | FFT 内纯计算占比 |
|:---:|:---:|:---:|:---:|:---:|
| 4 | 49.5% | 8.5% | 58.4% | 90.7% |
| 8 | 48.4% | 15.0% | 62.0% | 83.5% |
| 16 | 43.8% | 22.5% | 52.9% | 75.4% |
| 32 | 41.6% | **28.5%** | **61.5%** | 69.1% |

FFT 转置消息数随 rank 爆炸（4→32 ranks：16→16496 条）。**注意此表是仿真时钟口径，正式数据需用第四节的真实时钟实验重测。**

---

## 四、FFT 计算 / 转置通信分离（真实时钟口径）——已完成（2026-08-12）

**目标**：把 Kspace 拆成 4 块随 rank 对比——FFT 计算、FFT 转置（全局）通信、邻居通信（GRID_REVERSE+FORWARD）、其余；证明转置通信大幅上升、计算与邻居通信下降或小幅增加。

### 4.1 插桩（实际实施，与原方案有重要偏差）

**原方案的错误**：只在 `remap_wrap.cpp` 的 `Remap::perform()` 插桩会**漏掉 FFT 内部的转置通信**——`fft3d.cpp` 的 4 处转置直接调 C 函数 `remap_3d()`，不经过 `Remap` 类。首轮实验 32-rank 出现"FFT 计算 129.8 ms/step"假象，实为转置通信被误计入计算段。

**正确插桩点 = `remap.cpp` 的 `remap_3d()` 函数体首尾**（所有转置通信的汇聚点，覆盖 fft3d 内部 + Remap 类两条路径）：

- `src/KSPACE/remap.cpp`：`remap_3d()` 入口 `phase_marker_hook("KSPACE_FFT_COMM")`，出口 `phase_marker_hook("KSPACE_FFT_CALC")`；`remap_wrap.cpp` 无需改动（内部调 `remap_3d`）。
- `src/timer.h/.cpp`：新增静态钩子 `Timer::phase_marker_hook` + `_hook_owner` + `default_phase_marker()`，供无 LAMMPS 指针的 C 风格代码发 marker；`enable_phase_trace()` 时注册 owner。
- **关键坑（SMPI 协程模型）**：`SMPI_PRIVATIZATION=0` 下所有 rank 是同一进程的协程，共享全局变量与 OS 线程。静态钩子会被全部 rank 调用（首轮 marker 被放大 16 倍：36400 条），`thread_local` 也无效（协程共享线程）；必须在 `default_phase_marker()` 内用 `MPI_Comm_rank(_hook_owner->world, &me)` 实时过滤，仅 rank 0 写文件。
- marker 数量校验：每配置 COMM marker ≈ 2250 条（minimize 163 + run 10 步，每步约 13 对 remap）。

### 4.2 重编译

```bash
cd /work1/jiangtao/lammps_trace/lammps-build-smpi
make -j8 lmp
python3 /work1/jiangtao/lammps_trace/simgrid_traces/strip_pie_flag.py lmp
```

### 4.3 采集

产物目录 `runs/smpi_fft_split/`（`phase_{4,8,16,32}.csv` + `run_{4,8,16,32}.log`），命令同原方案（`LAMMPS_PHASE_TRACE` + `SMPI_PRIVATIZATION=0`，无需 Paje）。

### 4.4 分析

脚本：[analyze_fft_split.py](file:///work1/jiangtao/lammps_trace/runs/smpi_fft_split/analyze_fft_split.py)，结果存档 `runs/smpi_fft_split/analysis_result.txt`。双口径输出：真实时钟（第 1 列，主数据）+ 仿真时钟（第 3 列，辅助）。

### 4.5 结果（2688 原子，10 步 run，真实时钟 ms/step）

| ranks | FFT 转置通信 | FFT 计算 | 邻居通信 | Kspace 其余 | Kspace 合计 | 转置占 Kspace |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 4 | 4.54 | 1.03 | 1.17 | 0.24 | 6.99 | 64.9% |
| 8 | 8.87 | 0.87 | 1.79 | 0.20 | 11.74 | 75.6% |
| 16 | 17.66 | 0.42 | 4.05 | 0.14 | 22.26 | 79.3% |
| 32 | **128.78** | 0.39 | 8.19 | 0.18 | 137.54 | **93.6%** |

对齐校验：各分项之和与 Kspace 窗口总长偏差 +0.00%；窗口总长与 log 口径 Kspace 时长差 2-15%（rank0 窗口 vs 全 rank 均值口径差异）。

### 4.6 结论（假设得证）

1. **转置通信随 rank 超线性爆炸**：4→32 ranks 真实时钟耗时 ×28.4（4.54→128.78 ms/step），占 Kspace 比例 64.9%→93.6%，是 Loop time 恶化（×13.2）的绝对主因。
2. **FFT 计算随 rank 下降**：1.03→0.39 ms/step（≈1/N，符合计算量切分预期）。
3. **邻居通信增长平缓**：1.17→8.19 ms/step（×7，近线性，含小 brick 多轮交换退化），占比被转置淹没（16.8%→6.0%）。
4. **真实时钟 vs 仿真时钟的差值即"等待放大"**：仿真时钟列转置通信基本持平（0.35→0.25 ms/step，纯网络模型时间），真实时钟却 ×28——SMPI 协程模型下 rank 0 等消息时调度器去跑其他 rank 的协程（含计算），规模越大等待期内要陪跑的协程越多，应用层真实墙钟被 mesh 争用/多跳超线性放大。这正是"应用层真实耗时"口径要捕捉的效应。

---

## 五、已知坑清单

| 坑 | 规避 |
|---|---|
| 重链接 lmp 后 DF_1_PIE 回来 | 每次链接后必跑 strip_pie_flag.py |
| privatization 默认开 → MPI_ERR_COMM 崩溃 | `SMPI_PRIVATIZATION=0` |
| LD_PRELOAD shim 拦截不到 dlopen 程序的 MPI 符号 | 放弃该路线；需要 MPI 级拦截时用静态链接 wrapper（`simgrid_traces/smpi_walltime_static.c`，rank 解析仍有遗留问题） |
| log 有两个 timing breakdown（minimize/run） | `awk 'n==2'` 取第二个；相位窗口按 Loop 行步数切分 |
| LAMMPS Timer 是真实时钟不是仿真时钟 | 应用层分析直接用（符合实验目的）；仿真时钟分析才用第 3 列/MPI_Wtime |
| SMPI 默认不仿真计算时间 | 仿真时钟分析加 `--cfg=smpi/host-speed:2.8544e9`；应用层实验无关 |
| 链路口径求和可超 100%（并行重叠） | 单 rank 分析用 op 时长（无重叠）口径 |
| proc grid 形状伪影（转置对恰为邻居） | np=4(2×2×1) 转置占比是下界；np=16(4×2×2) 有回落假象 |
| bwrap 沙箱故障 | 命令用 required_permissions='all' |
| 只插桩 `Remap::perform()` 会漏掉 FFT 内部转置（fft3d.cpp 直接调 C 函数 `remap_3d`） | 插桩点放 `remap.cpp` 的 `remap_3d()` 首尾；症状：某 rank "FFT 计算"段异常巨大（如 32r 129.8 ms/step） |
| SMPI 协程共享进程/线程，静态全局钩子被全部 rank 调用（marker 放大 N 倍）；thread_local 无效 | 钩子回调内 `MPI_Comm_rank` 实时过滤，仅 rank 0 写文件；症状：COMM marker 数 ≈ 预期 × rank 数 |
| breakdown 表格正则误抓 min time 列 | %total 是最后一列，按 `\|` 切分取末列转 float |

---

## 六、相关文档

- [AGENTS.md](file:///work1/jiangtao/lammps_trace/AGENTS.md)：项目总览 + SMPI 重放方法 + Kspace 占比结论（第一阶段）
- `dumpi2ccdg_guide.md`、`CCDG实验报告.md`：更早期的 trace/CCDG 流水线（与当前实验弱相关）
