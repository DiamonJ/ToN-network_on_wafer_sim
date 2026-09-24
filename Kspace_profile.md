# Kspace 操作细节耗时提取方法（Kspace Profile）

> 本文档总结在 SMPI 2D mesh 虚拟平台上，如何从 LAMMPS 应用层真实耗时中提取 Kspace 相位内部的操作级细节（FFT 计算 / FFT 转置通信 / 邻居通信），包含插桩机制、运行方法、分析流程与已知注意事项。

---

## 一、实验目的与数据口径

### 1.1 研究问题与假设

研究**强可扩展性下 2D mesh 的通信特征**：随着 rank 增加，Kspace 为何成为主要开销。

- **假设**：2D mesh 随规模增加通信性能下降，Kspace 中需要全局通信的 FFT（转置通信）成为通信开销的主要来源。
- **待证明**：Kspace 内部，**计算和邻居通信随 rank 增加而降低或小幅增加，而 FFT 全局（转置）通信随 rank 大幅增加**。
- **判据**：转置通信 per-step 随 rank 大幅上升，FFT 计算随 rank 下降（≈1/N），邻居通信持平或小幅变化 → 假设得证。

### 1.2 拆分目标

把每步的 Kspace 相位耗时拆成 4 块随 rank 数对比：

| 分项 | 含义 |
|---|---|
| FFT 转置通信 | brick 分解 ↔ FFT 分解的数据重排通信（全局，强扩展瓶颈） |
| FFT 计算 | FFTW 的 1D FFT 变换及 brick2fft 准备代码 |
| 邻居通信 | GRID_REVERSE（ρ 密度幽灵交换）+ GRID_FORWARD（E 电场幽灵交换），仅 brick 6 面邻居 |
| Kspace 其余 | FORCE（场力插值）、REDUCE（Allreduce）及窗口残余 |

### 1.3 数据口径（用户要求）

LAMMPS 作为应用真实运行，rank 间通信由 SimGrid 按 2D mesh 模型计时，耗时取 **LAMMPS 应用层真实时钟**（`platform::walltime()` = `std::chrono::steady_clock`，不是 `MPI_Wtime` 仿真时钟）。SMPI 下进程按仿真网络结果真实阻塞，因此 breakdown 即"应用层真实耗时"。

### 1.4 问题发现过程（三阶段递进）

**阶段一：发现 Kspace 占比随 rank 单调上升（相位级粗拆）**

SMPI 2D mesh 重放真实 LAMMPS（LiAlOCl 算例，lj/cut/coul/long + PPPM），log 的 `MPI task timing breakdown` 显示：

| ranks | Loop time (s) | Kspace %total | Comm %total | Pair %total |
|:---:|:---:|:---:|:---:|:---:|
| 4 | 0.1495 | 72.5% | 18.7% | 6.2% |
| 8 | 0.1596 | 78.3% | 17.5% | 2.3% |
| 16 | 0.2913 | 80.7% | 16.9% | 0.7% |
| 32 | 1.5252 | **92.9%** | 6.3% | 0.1% |

Loop time 超线性恶化（rank×8 → 时间×10.2），增量几乎全来自 Kspace。**但此时只能定位到"Kspace 是祸首"，尚不知道 Kspace 内部是计算还是哪类通信主导**——需要更细的拆分。

**阶段二：仿真时钟口径初步定位 FFT 转置（Paje 二级解剖）**

用 Paje trace（仿真时钟 + host-speed 口径）把 Kspace 内的通信事件归位到子相位，发现：

| ranks | FFT 阶段占 run 段 | 转置通信占 FFT 内 | FFT 内纯计算占比 |
|:---:|:---:|:---:|:---:|
| 4 | 49.5% | 8.5% | 90.7% |
| 8 | 48.4% | 15.0% | 83.5% |
| 16 | 43.8% | 22.5% | 75.4% |
| 32 | 41.6% | **28.5%** | 69.1% |

FFT 转置消息数随 rank 爆炸（4→32 ranks：16→16496 条），方向符合假设；**但此表是仿真时钟口径，不是用户要求的应用层真实耗时，只能作佐证**——需要在真实时钟口径下重新分离 FFT 计算与转置通信，即本文档的方法。

**阶段三：真实时钟口径细拆（本文档），两次插桩异常的发现与修复**

首轮实验（按原方案只插桩 `Remap::perform()`）即暴露异常：

1. **异常 A：32-rank "FFT 计算"达 129.8 ms/step**（预期应 ≈1/N 下降），且每步仅 1 对 COMM/CALC marker（预期 ~13 对）。逐段检查发现 124.7 ms 是一个连续的 CALC 段 → 追查调用链：`fft3d.cpp` 的 4 处转置直接调 C 函数 `remap_3d()`，不经过 `Remap` 类，转置通信时间被误计进了计算段。**修复：插桩点改到 `remap.cpp` 的 `remap_3d()` 首尾**（所有转置路径的汇聚点，见第三节）。
2. **异常 B：COMM marker 暴增到 36400 条**（预期 2275 条，恰为 16 倍）。排查：`SMPI_PRIVATIZATION=0` 下所有 rank 是同一进程的协程，静态钩子被全部 rank 调用；先试 `thread_local` 无效（协程共享 OS 线程），**最终用钩子内 `MPI_Comm_rank` 实时过滤解决**（见第三节 3.4）。
3. 并行发现 log 解析正则误抓 breakdown 表格的 min time 列（%total 实为最后一列）。

修复后数据闭合（拆分和偏差 +0.00%，marker 数恢复预期），得到第六节的正式结果，假设得证。

---

## 二、原理：源码级相位 marker 机制

### 2.1 为什么用 marker 而不是消息特征

- DUMPI/消息层面，kspace 内所有通信都用同一个 `world` 通信域，comm 字段无法区分局部/全局；
- 消息大小阈值法在强扩展下失效（所有消息随 rank 数变小）；
- **marker 按代码位置归类**："这段代码在做什么"无歧义，不依赖任何消息特征。

### 2.2 marker 写入机制

- `src/timer.cpp` 支持环境变量 `LAMMPS_PHASE_TRACE=<csv路径>`：仅 rank 0 将相位 marker 写入 CSV。
- CSV 格式（每行 3 列）：

```
walltime_ns,phase_name,sim_time_s
```

| 列 | 时钟 | 用途 |
|---|---|---|
| 第 1 列 `walltime_ns` | CLOCK_MONOTONIC 真实时钟（ns） | **主口径**，应用层真实耗时分析 |
| 第 3 列 `sim_time_s` | `MPI_Wtime`（SMPI 下为仿真时钟） | 辅助口径，剔除"对端计算拖慢"后的纯通信/自身时间 |

### 2.3 两类 marker 语义（重建时间线的关键）

| 类型 | 行名 | 语义 |
|---|---|---|
| 主相位 | `Pair` / `Kspace` / `Comm` / … | **END marker**：`(上一主相位行 ts, 本行 ts]` 属于该相位 |
| 子相位 | `KSPACE_*` | **START marker**：段持续到下一个任意 marker |

### 2.4 每步 Kspace 内的 marker 序列（16-rank 实测）

```
GRID_REVERSE → FFT → (COMM → CALC) × 13 → GRID_FORWARD → FORCE → REDUCE
```

每对 `(COMM, CALC)` 对应一次 `remap_3d()` 调用（转置通信 + 其后 FFT 计算），每步 13 次 remap（正反向 FFT 各含 fft3d 内 4 处 + brick2fft 路径）。

---

## 三、插桩点清单

### 3.1 主相位 marker（LAMMPS 原生 timer，无需修改）

`Timer::stamp()` 在各相位结束时写入主相位行，自动产生 Kspace 窗口。

### 3.2 Kspace 子相位 marker（`src/KSPACE/pppm.cpp` 的 `compute()`）

| marker | 埋点位置 | 语义 |
|---|---|---|
| `KSPACE_GRID_REVERSE` | ρ 网格逆向交换前（L639 附近） | 邻居通信 |
| `KSPACE_FFT` | brick2fft + poisson + fft2brick 前（L642） | FFT 阶段开始 |
| `KSPACE_GRID_FORWARD` | E 场正向交换前（L655） | 邻居通信 |
| `KSPACE_FORCE` | 场力插值前（L676） | 其余 |
| `KSPACE_REDUCE` | MPI_Allreduce 前（L685） | 其余 |

均通过 `timer->mark_subphase(name)` 写入。

### 3.3 FFT 转置 marker（`src/KSPACE/remap.cpp` 的 `remap_3d()` 首尾）

```cpp
void remap_3d(FFT_SCALAR *in, FFT_SCALAR *out, FFT_SCALAR *buf,
              struct remap_plan_3d *plan)
{
  LAMMPS_NS::Timer::phase_marker_hook("KSPACE_FFT_COMM");   // 入口：转置通信段
  ...（p2p 或 Alltoallv 两个分支的通信代码）...
  LAMMPS_NS::Timer::phase_marker_hook("KSPACE_FFT_CALC");   // 出口：其后为 FFT 计算
}
```

**为什么必须插在 `remap_3d()` 而不是 `Remap::perform()`**：转置通信有两条路径——

1. `Remap` 类（`remap_wrap.cpp`）：pppm.cpp 的 `remap->perform()`（brick2fft）；
2. C 函数 `remap_3d()`：`fft3d.cpp` 内 4 处直接调用（pre/mid1/mid2/post remap）。

`Remap::perform()` 内部也调 `remap_3d()`，所以 **`remap_3d()` 是所有转置通信的唯一汇聚点**。只插桩 `Remap::perform()` 会漏掉 FFT 内部的转置（症状：某 rank 的"FFT 计算"段异常巨大，如 32-rank 首轮出现 129.8 ms/step 假象）。

### 3.4 钩子机制（`src/timer.h/.cpp`）

`remap.cpp`/`fft3d.cpp` 是 C 风格代码，无 LAMMPS 指针，通过静态钩子发 marker：

```cpp
// timer.h（Timer 类内）
static void (*phase_marker_hook)(const char *);
static void default_phase_marker(const char *);
static class Timer *_hook_owner;   // enable_phase_trace() 时注册为 rank 0 的 Timer

// timer.cpp
void Timer::default_phase_marker(const char *name)
{
  if (!_hook_owner) return;
  int me = 0;
  MPI_Comm_rank(_hook_owner->world, &me);   // 关键：SMPI 协程模型下必须实时过滤
  if (me == 0) _hook_owner->mark_subphase(name);
}
```

**SMPI 协程陷阱**：`SMPI_PRIVATIZATION=0` 下所有 rank 是同一进程的协程，共享全局变量与 OS 线程。静态钩子会被全部 rank 调用（marker 被放大 N 倍，实测 16-rank 时 36400 条 vs 预期 2275 条），`thread_local` 也无效（协程共享线程）。必须在回调内用 `MPI_Comm_rank` 实时判断，仅 rank 0 写文件。

---

## 四、运行方法

### 4.1 前置条件（缺一不可）

| 条件 | 原因 |
|---|---|
| `SMPI_PRIVATIZATION=0` | privatization 会使 `MPI_COMM_WORLD` 绑定为 NULL，报 `MPI_ERR_COMM` |
| 二进制清 `DF_1_PIE` 标志 | glibc 拒绝 dlopen 带 PIE 标志的可执行文件；每次重链接后必跑 `python3 simgrid_traces/strip_pie_flag.py lammps-build-smpi/lmp` |
| `-Wl,--export-dynamic-symbol=main` | SMPI dlopen 后需从动态符号表找 `main`（现有构建已满足） |
| `timer normal`（默认）以上 | `timer loop/off` 会静默丢失相位标记 |

### 4.2 重编译（修改插桩后）

```bash
cd /work1/jiangtao/lammps_trace/lammps-build-smpi
make -j8 lmp
python3 /work1/jiangtao/lammps_trace/simgrid_traces/strip_pie_flag.py lmp
```

### 4.3 采集（4/8/16/32 ranks，2688 原子算例）

**推荐：封装脚本 [run_kspace_profile.sh](file:///work1/jiangtao/lammps_trace/run_kspace_profile.sh)**（自动完成 PIE 检查、mesh 尺寸选择/平台生成、输入文件拷贝、运行、分析归档）：

```bash
./run_kspace_profile.sh \
    -i cases/lialocl_coul_2688/in.lammps \
    -d cases/lialocl_coul_2688/data.LiAlOCl_nvt_charge \
    -n 16            # rank 数；可选 -s 步数（默认 10）-t 超时秒（默认 1500）
```

产物目录 `runs/<时间戳>_<np>ranks/`：过程文件（`phase_<np>.csv`、`run_<np>.log`、`log.lammps`）在根目录，最终结果在 `result/`（`log.lammps` + `kspace_breakdown.txt`）。mesh 尺寸自动按 np 整除分解选择（4→2×2、8→2×4、16→4×4、32→4×8），无现成平台 XML 时在运行目录自动生成。

**batch 模式**（多算例 × 多 rank 批量，`-b` 指定 batch 文件）：

```bash
./run_kspace_profile.sh -b batch.txt        # 可选 -s 步数 -t 超时秒
```

batch 文件格式（每行一个算例，`#` 开头为注释）：

```
<in.lammps 路径>,<data 文件路径>,[r1,r2,...][,步数]
cases/lialocl_coul_2688/in.lammps,cases/lialocl_coul_2688/data.LiAlOCl_nvt_charge,[4,8,16,32]
cases/lialocl_coul_1344/in.lammps,cases/lialocl_coul_1344/data.LiAlOCl_nvt_charge,[4,8],20
```

行内步数可选，缺省用全局 `-s`（默认 10）。

产物目录 `runs/<时间戳>_batch/`：每个 (算例, rank) 的过程文件（含记录 mesh 结构的 `mesh.txt`）与 `result/{log.lammps,kspace_breakdown.txt}` 在 `<算例名>/<np>ranks/` 下；根目录的 **`batch_result.csv`** 汇总所有 rank 实验的 mesh 结构、MPI breakdown（7 相位的 avg 时长 + %total）与 Kspace breakdown（4 分项 ms/step + 占 Kspace 百分比）。单个 rank 失败不中断批次，在 CSV 的 status 列标记。

等价的裸命令（多配置批量）：

```bash
mkdir -p /work1/jiangtao/lammps_trace/runs/smpi_fft_split
cd /work1/jiangtao/lammps_trace/cases/lialocl_coul_2688
for cfg in "platform_mesh_2x2.xml 4" "platform_mesh_2x4.xml 8" "platform_mesh_4x4.xml 16" "platform_mesh_4x8.xml 32"; do
  set -- $cfg
  SMPI_PRIVATIZATION=0 \
  LAMMPS_PHASE_TRACE=/work1/jiangtao/lammps_trace/runs/smpi_fft_split/phase_$2.csv \
    timeout 1500 smpirun \
    -platform /work1/jiangtao/lammps_trace/simgrid_traces/$1 \
    -np $2 /work1/jiangtao/lammps_trace/lammps-build-smpi/lmp \
    -in in.lammps -v STEPS 10 \
    > /work1/jiangtao/lammps_trace/runs/smpi_fft_split/run_$2.log 2>&1
  echo "np=$2 rc=$?"
done
```

**无需 Paje tracing**，纯应用层运行。

**marker 数量校验**：每配置 `grep -c KSPACE_FFT_COMM phase_*.csv` ≈ 2250 条（minimize 163 + run 10 步，每步约 13 对 remap）。若 ≈ 预期 × rank 数，说明钩子 rank 过滤失效。

---

## 五、分析方法

**脚本**：[runs/smpi_fft_split/analyze_fft_split.py](file:///work1/jiangtao/lammps_trace/runs/smpi_fft_split/analyze_fft_split.py)

```bash
python3 runs/smpi_fft_split/analyze_fft_split.py runs/smpi_fft_split
# 可选：指定目录与 ranks
python3 analyze_fft_split.py <结果目录> 4,8,16,32
```

### 5.1 处理流程

1. **切 run 阶段**：从 `run_*.log` 最后一个 `Loop time of X on N procs for M steps` 取 run 步数 M（算例含 minimize 100 步 + run 10 步，log 里有两个 breakdown，取第二个）；Kspace 窗口序列的**最后 M 个**属于 run 阶段。
2. **重建子段**：run 区间内每个 `KSPACE_*` marker 为一段起点，段终点 = 下一个任意 marker（主相位行自然结束子段），按 Kspace 窗口裁切。
3. **分类累加**（BUCKET 映射）：

| 段名 | 归类 |
|---|---|
| `KSPACE_FFT_COMM` | FFT 转置通信 |
| `KSPACE_FFT_CALC`、`KSPACE_FFT`（含残余） | FFT 计算 |
| `KSPACE_GRID_REVERSE`、`KSPACE_GRID_FORWARD` | 邻居通信 |
| `KSPACE_FORCE`、`KSPACE_REDUCE`、窗口未覆盖残余 | Kspace 其余 |

4. **双时钟口径输出**：第 1 列（真实时钟，主数据）+ 第 3 列（仿真时钟，辅助），每步归一为 ms/step。
5. **对齐校验**：
   - 各分项之和 vs Kspace 窗口总长（偏差应为 ~0%，否则有暗时间或 marker 丢失）；
   - 窗口总长 vs log 第二个 breakdown 的 `Kspace %total × Loop time`（差 2-15% 属正常，rank 0 窗口 vs 全 rank 均值口径差异；注意 %total 是表格**最后一列**，按 `|` 切分取末列，正则抓首列会误取 min time）。

### 5.2 两个时钟口径的解读

| 口径 | 含义 | 特征 |
|---|---|---|
| 真实时钟（第 1 列，**主数据**） | rank 0 感受到的应用层墙钟 | 转置通信 4→32 ranks ×28.4 |
| 仿真时钟（第 3 列，辅助） | SMPI 网络模型时间（默认不计计算） | 转置通信基本持平（0.35→0.25 ms/step） |

两者差值 = **"等待放大"**：SMPI 协程模型下 rank 0 等消息时，调度器去跑其他 rank 的协程（含计算）；规模越大，等待期内要陪跑的协程越多，应用层真实墙钟被 mesh 争用/多跳超线性放大。这正是"应用层真实耗时"口径要捕捉的效应。

---

## 六、实测结果（2688 原子，10 步 run，真实时钟 ms/step）

| ranks | FFT 转置通信 | FFT 计算 | 邻居通信 | Kspace 其余 | Kspace 合计 | 转置占 Kspace |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 4 | 4.54 | 1.03 | 1.17 | 0.24 | 6.99 | 64.9% |
| 8 | 8.87 | 0.87 | 1.79 | 0.20 | 11.74 | 75.6% |
| 16 | 17.66 | 0.42 | 4.05 | 0.14 | 22.26 | 79.3% |
| 32 | **128.78** | 0.39 | 8.19 | 0.18 | 137.54 | **93.6%** |

**结论**：
1. **转置通信超线性爆炸**（×28.4），占 Kspace 比例 64.9%→93.6%，是 Loop time 恶化（×13.2）的绝对主因；
2. **FFT 计算随 rank 下降** ≈1/N（1.03→0.39 ms/step）；
3. **邻居通信近线性缓增**（×7，含小 brick 多轮交换退化），占比被转置淹没（16.8%→6.0%）。

---

## 七、已知坑清单

| 坑 | 症状 | 规避 |
|---|---|---|
| 只插桩 `Remap::perform()` 漏掉 FFT 内部转置（fft3d.cpp 直调 C 函数 `remap_3d`） | 某 rank "FFT 计算"段异常巨大（32r 129.8 ms/step） | 插桩点放 `remap.cpp` 的 `remap_3d()` 首尾 |
| SMPI 协程共享进程/线程，静态全局钩子被全部 rank 调用；`thread_local` 无效 | COMM marker 数 ≈ 预期 × rank 数 | 钩子回调内 `MPI_Comm_rank` 实时过滤，仅 rank 0 写文件 |
| 重链接 lmp 后 `DF_1_PIE` 标志回来 | SMPI dlopen 失败 | 每次链接后必跑 `strip_pie_flag.py` |
| privatization 默认开 | `MPI_ERR_COMM` 崩溃 | `SMPI_PRIVATIZATION=0` |
| log 有两个 timing breakdown（minimize/run） | 取错阶段数据 | `awk '/MPI task timing breakdown/{n++} n==2'` 取第二个；相位窗口按 Loop 行步数切分 |
| breakdown 表格正则误抓 min time 列 | Kspace %total 解析成 ~0 | %total 是最后一列，按 `\|` 切分取末列转 float |
| `timer loop/off` | 相位 marker 静默丢失 | 输入脚本保持 `timer normal`（默认）以上 |
| LAMMPS Timer 是真实时钟不是仿真时钟 | — | 应用层分析直接用（符合本实验目的）；仿真时钟分析用第 3 列 |

---

## 八、相关文件索引

| 文件 | 说明 |
|---|---|
| `lammps-src/src/timer.h/.cpp` | marker 写入机制 + 静态钩子 |
| `lammps-src/src/KSPACE/pppm.cpp` | 5 个主 marker 埋点 |
| `lammps-src/src/KSPACE/remap.cpp` | 转置 COMM/CALC marker（`remap_3d()` 首尾） |
| `runs/smpi_fft_split/analyze_fft_split.py` | 双口径分析脚本 |
| `runs/smpi_fft_split/phase_{4,8,16,32}.csv` | 相位 trace 原始数据 |
| `runs/smpi_fft_split/run_{4,8,16,32}.log` | 运行日志（含 breakdown） |
| `runs/smpi_fft_split/analysis_result.txt` | 分析结果存档 |
| `simgrid_traces/platform_mesh_{2x2,2x4,4x4,4x8}.xml` | 2D mesh 平台（host 1.42 Gflops，链路 6.8 GB/s，延迟 10 ns） |
| `Kspace_FFT实验手册.md` | 实验总手册（环境、历史结果、本文档的上游文档） |
