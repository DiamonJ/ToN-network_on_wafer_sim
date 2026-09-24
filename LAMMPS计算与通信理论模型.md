# LAMMPS 计算与通信流程及理论模型

## 1. 文档范围

本文面向 `run_noc_pipeline.sh` 生成的 LAMMPS 输入，说明：

1. LAMMPS 按什么顺序执行计算；
2. 每个计算阶段的理论工作量如何估算；
3. LAMMPS 在什么位置发生通信；
4. 每个通信阶段的消息量和字节数如何估算；
5. 如何得到每个 rank 的通信量 T1 和计算量 C1。

当前覆盖的输入场景：

```text
Cu：EAM
H2O：LJ + Coulomb
LiAlOCl：LJ + Coulomb
short：截断库仑
long：截断实空间 + PPPM Kspace/FFT
processors：K × K × 1
boundary：p p p
integrator：NVE/Velocity Verlet
```

源码入口：

```text
lammps-src/src/verlet.cpp
lammps-src/src/comm_brick.cpp
lammps-src/src/neighbor.cpp
lammps-src/src/pair_lj_cut_coul_cut.cpp
lammps-src/src/MANYBODY/pair_eam.cpp
lammps-src/src/KSPACE/pppm.cpp
lammps-src/src/grid3d.cpp
lammps-src/src/KSPACE/remap.cpp
```

## 2. T1 和 C1 的定义

本文定义：

```text
T1(r) = 一个正式 timestep 中 rank r 发出的应用层 payload bytes

C1(r) = 一个正式 timestep 中 rank r 执行的理论计算工作量
```

推荐同时输出以下 C1 指标：

```text
本地原子数
neighbor list entries
真正进入 force cutoff 的 pair 数
PPPM 网格更新数
FFT 工作量
模型 scalar operations
```

原因是“理论 FLOPs”和“硬件运行周期”不是同一个量：

- 除法、平方根与加法的硬件代价不同；
- 编译器可能使用 SIMD 和 FMA；
- cache miss 和访存时间不计入 FLOPs；
- DUMPI 中由时间换算出的 `compute_ops` 不是真实硬件 FLOPs。

此外要区分：

```text
T1_send：rank 发出的唯一 payload
T1_recv：rank 接收的 payload
T1_nic：send + recv
T1_hop：Σ(message_bytes × mesh_hops)
```

本文默认使用 `T1_send`，因为它与当前 `wse_plan.json` 和 CCDG SEND 节点口径一致。

## 3. 输入中决定计算量和通信量的参数

### 3.1 几何与原子

```text
boundary
processors
lattice / region / create_atoms
read_data
replicate
atom_style
原子坐标、类型、质量、电荷
```

它们决定：

- 全局盒子尺寸；
- 原子密度；
- 每个 rank 的 owned atoms；
- 每个 rank 的 ghost atoms；
- pair 和 neighbor list 的实际数量。

### 3.2 势函数

```text
pair_style
pair_coeff
kspace_style
kspace accuracy
```

它们决定：

- force cutoff；
- 每个 pair 的计算 kernel；
- 是否执行 PPPM；
- PPPM 网格规模和 FFT 工作量；
- EAM 是否产生额外 pair forward/reverse 通信。

### 3.3 邻居表

```text
neighbor <skin> bin
neigh_modify delay <D> every <E> check <yes|no>
```

通信 ghost 宽度通常为：

```math
g = r_{\mathrm{cut,max}} + s_{\mathrm{skin}}
```

邻居表 cutoff 也是：

```math
r_{\mathrm{neigh}} = r_{\mathrm{cut,max}} + s_{\mathrm{skin}}
```

### 3.4 时间积分和输出

```text
timestep
velocity
fix nve
thermo
run
```

它们决定：

- 积分次数；
- 原子移动速度；
- neighbor rebuild 和 atom migration 的可能频率；
- 哪些 timestep 需要能量、virial 和 collective。

## 4. 总体运行顺序

LAMMPS 分为一次 setup 和多个正式 timestep。

```text
Setup
  → 建立处理器拓扑
  → exchange
  → borders
  → neighbor build
  → 初始 pair/Kspace force
  → 初始 reverse communication

每个 timestep
  → initial_integrate
  → neighbor->decide
      ├→ 不重建：forward_comm
      └→ 重建：exchange + borders + neighbor build
  → pair compute
  → Kspace compute（仅 long）
  → reverse_comm
  → final_integrate
  → output（按 thermo/output 周期）
```

setup 和 run 必须分开统计。当前 pipeline 的 trim-only CCDG 和
`wse_plan scope="run"` 只保留正式 timestep。

## 5. Setup 阶段

## 5.1 处理器网格和子域

设：

```text
全局盒子：Lx × Ly × Lz
处理器网格：Px × Py × Pz
rank 数：R = Px × Py × Pz
```

均匀分解时，每个 rank 子域尺寸为：

```math
l_x = L_x/P_x,\quad
l_y = L_y/P_y,\quad
l_z = L_z/P_z
```

全局原子密度：

```math
\rho = N/(L_xL_yL_z)
```

平均 owned atoms：

```math
N_{\mathrm{local,avg}} = N/R
```

精确 `Nlocal(r)` 需要按坐标判断原子落在哪个 rank 子域。

当前 pipeline 使用：

```text
processors K K 1
```

因此：

```text
X、Y 维度跨 rank
Z 维度全部由本 rank 持有
Z 周期 ghost 可以本地复制，不产生 Z 网络通信
```

## 5.2 CommBrick setup

`CommBrick::setup()` 计算：

```text
ghost width
每个维度需要传播多少轮
每轮 send_proc / recv_proc
PBC 标记
发送 slab 范围
```

维度 d 的近似传播轮数：

```math
m_d = \left\lfloor g/l_d \right\rfloor + 1
```

总 swap 数：

```math
N_{\mathrm{swap}} = 2(m_x+m_y+m_z)
```

只有 `P_d>1` 且邻居不是本 rank 的 swap 才产生网络通信。

## 5.3 Initial exchange

`exchange()` 将位置已经不属于本 rank 子域的原子迁移给其他 rank。

若方向 d 有 `N_migrate(r,d)` 个原子：

```math
T_{\mathrm{exchange}}(r,d)
= N_{\mathrm{migrate}}(r,d)
\times S_{\mathrm{exchange}}\times 8
```

其中 `S_exchange` 是每个迁移原子打包的 double 数。

迁移量依赖实际坐标和速度，不能只从原子总数精确推导。

## 5.4 Borders

`borders()` 建立 ghost atom 并保存每个 swap 的 `sendlist`。

对于均匀密度，维度 d、第 m 轮的 slab 厚度可近似为：

```math
h_{d,m}=\max(0,\min(l_d,g-ml_d))
```

在依次完成前面维度通信后，有效截面积：

```math
A_d=
\prod_{j<d}\min(L_j,l_j+2g)
\times
\prod_{j>d}l_j
```

单方向发送原子数：

```math
N_{\mathrm{send}}(r,d,m)
\approx \rho h_{d,m}A_d
```

两个方向的总发送 atom-copies：

```math
N_{\mathrm{copies}}(r)
\approx
\sum_{d:P_d>1}\sum_m 2\rho h_{d,m}A_d
```

常见 border pack 大小：

```text
atom_style atomic：6 doubles/atom
atom_style charge：7 doubles/atom
```

所以：

```math
T_{\mathrm{borders}}(r)
=N_{\mathrm{copies}}(r)
\times S_{\mathrm{border}}\times 8
```

精确计算时不能只使用连续密度公式，必须像 LAMMPS 一样按 swap 顺序对真实坐标
进行 slab 筛选，并将前一维收到的 ghost 加入后一维候选。

## 5.5 Neighbor build

neighbor build 的主要计算：

```text
将 owned+ghost atoms 放入空间 bins
遍历每个 owned atom 周围的 stencil bins
计算候选原子距离
将满足 neighbor cutoff 的原子写入 neighbor list
```

均匀体系中，单原子的理论邻居数：

```math
n_{\mathrm{neigh}}
\approx \rho\frac{4\pi}{3}(r_{\mathrm{cut}}+s_{\mathrm{skin}})^3
```

half neighbor list 且 Newton pair 开启时，每 rank 的 list entries：

```math
N_{\mathrm{list}}(r)
\approx
\frac{1}{2}N_{\mathrm{local}}(r)n_{\mathrm{neigh}}
```

neighbor build 工作量：

```math
C_{\mathrm{neighbor}}(r)
\approx
F_{\mathrm{bin}}(N_{\mathrm{local}}+N_{\mathrm{ghost}})
+F_{\mathrm{candidate}}N_{\mathrm{candidate}}
+F_{\mathrm{write}}N_{\mathrm{list}}
```

## 6. 正式 timestep：积分

## 6.1 Initial integrate

Velocity Verlet 的第一阶段近似为：

```math
\mathbf v(t+\Delta t/2)
=\mathbf v(t)+\frac{\Delta t}{2m}\mathbf f(t)
```

```math
\mathbf x(t+\Delta t)
=\mathbf x(t)+\Delta t\mathbf v(t+\Delta t/2)
```

三个空间维度中，每个原子包含速度和位置的乘加操作。

## 6.2 Final integrate

力计算完成后：

```math
\mathbf v(t+\Delta t)
=\mathbf v(t+\Delta t/2)
+\frac{\Delta t}{2m}\mathbf f(t+\Delta t)
```

可以使用近似：

```math
C_{\mathrm{integrate}}(r)
\approx 18N_{\mathrm{local}}(r)
```

这里将一次乘法和一次加法分别计为一个 scalar operation。

## 7. 正式 timestep：邻居判断和通信分支

`neighbor->decide()` 决定当前步走哪条路径。

## 7.1 不重建邻居表

```text
forward_comm
→ pair/Kspace
→ reverse_comm
```

这是稳态 timestep。

## 7.2 重建邻居表

```text
exchange
→ borders
→ neighbor build
→ pair/Kspace
→ reverse_comm
```

这是 rebuild timestep。

平均每步成本可以摊销为：

```math
T_{\mathrm{avg}}
=T_{\mathrm{steady}}
+T_{\mathrm{rebuild,extra}}/I_{\mathrm{rebuild}}
```

```math
C_{\mathrm{avg}}
=C_{\mathrm{steady}}
+C_{\mathrm{rebuild,extra}}/I_{\mathrm{rebuild}}
```

其中 `I_rebuild` 是平均多少步重建一次 neighbor list。

## 8. Forward communication

`forward_comm()` 将 owner atom 的最新坐标发送到 ghost 副本。

当前场景通常发送：

```text
x、y、z：3 doubles/atom
```

第 i 个 swap：

```math
T_{\mathrm{forward}}(r,i)
=N_{\mathrm{send}}(r,i)\times3\times8
```

每 rank：

```math
T_{\mathrm{forward}}(r)
=\sum_i N_{\mathrm{send}}(r,i)\times24
```

`forward_comm()` 复用最近一次 `borders()` 建立的 sendlist，不会每个 timestep
重新选择 ghost。

## 9. Pair force 计算

## 9.1 通用 pair 数

真实进入力计算 cutoff 的平均邻居数：

```math
n_{\mathrm{force}}
\approx\rho\frac{4\pi}{3}r_{\mathrm{cut}}^3
```

half neighbor list：

```math
N_{\mathrm{pair}}(r)
\approx\frac{1}{2}
N_{\mathrm{local}}(r)n_{\mathrm{force}}
```

多类型体系中，应按类型比例和 pair cutoff 求和：

```math
N_{\mathrm{pair}}(r)
\approx
\frac{1}{2}N_{\mathrm{local}}(r)
\sum_a\sum_b x_ax_b
\rho\frac{4\pi}{3}r_{c,ab}^3
```

## 9.2 LJ + Coulomb short

每个 neighbor candidate 首先执行：

```text
delx、dely、delz
rsq = delx² + dely² + delz²
cutoff 判断
```

可以近似为：

```text
Fcandidate ≈ 8 scalar operations
```

落入 cutoff 后执行：

```text
r2inv
Coulomb force
LJ r6 和 force
合并 fpair
更新 i、j 的三个力分量
可选 energy/virial tally
```

参数化模型：

```math
C_{\mathrm{pair}}(r)
=F_{\mathrm{candidate}}N_{\mathrm{list}}(r)
+F_{\mathrm{kernel}}N_{\mathrm{pair}}(r)
```

当前 scalar-equivalent 初值可以取：

```text
Fcandidate = 8
Fkernel = 45～50
```

`Fkernel` 必须通过 profiler 校准，不应当成硬件精确 FLOPs。

## 9.3 LJ + Coulomb long

实空间 pair 部分仍使用上述模型，但 Coulomb 实空间 kernel 与 short 略有不同。

额外的长程静电由 PPPM 执行：

```math
C_{\mathrm{long}}=
C_{\mathrm{real-space}}+C_{\mathrm{PPPM}}
```

## 9.4 EAM

EAM 计算顺序：

```text
遍历 neighbor pairs，累计电子密度
→ pair_reverse
→ 每个原子计算 embedding
→ pair_forward
→ 再次遍历 neighbor pairs，计算力
```

计算量：

```math
C_{\mathrm{EAM}}(r)
=F_{\rho}N_{\mathrm{pair}}(r)
+F_{\mathrm{embed}}
(N_{\mathrm{local}}+N_{\mathrm{ghost}})
+F_{\mathrm{force}}N_{\mathrm{pair}}(r)
```

EAM 的两个附加通信通常每 ghost atom 发送一个 double：

```math
T_{\mathrm{pair-forward}}(r)
\approx 8\sum_iN_{\mathrm{send}}(r,i)
```

```math
T_{\mathrm{pair-reverse}}(r)
\approx 8\sum_iN_{\mathrm{reverse}}(r,i)
```

## 10. PPPM 计算与通信

long 模式执行：

```text
particle_map
→ make_rho
→ grid_reverse
→ brick2fft
→ forward FFT
→ Poisson
→ inverse FFT
→ grid_forward
→ fieldforce
→ Allreduce
```

设：

```text
PPPM order：p
全局网格：Nx × Ny × Nz
全局网格点数：M = NxNyNz
```

## 10.1 电荷映射

每个原子影响约 `p³` 个网格点：

```math
C_{\mathrm{rho}}(r)
\approx F_{\mathrm{rho-grid}}
N_{\mathrm{local}}(r)p^3
```

## 10.2 Grid reverse

ghost grid 上的电荷密度需要累加回 owner。

```math
T_{\mathrm{grid-reverse}}(r)
=N_{\mathrm{reverse-grid}}(r)
\times1\times S_{\mathrm{FFT}}
```

double FFT precision 时：

```text
SFFT = 8 bytes
```

## 10.3 FFT remap

brick decomposition 与 FFT pencil decomposition 之间需要数据转置。

对 source rank s 和 destination rank d：

```math
V_{\mathrm{overlap}}(s,d)
=
V(\mathrm{source\ extent}_s
\cap
\mathrm{destination\ extent}_d)
```

消息大小：

```math
T_{\mathrm{remap}}(s,d)
=V_{\mathrm{overlap}}(s,d)
\times n_{\mathrm{qty}}\times S_{\mathrm{FFT}}
```

其中：

```text
实数 brick-to-FFT：nqty=1
复数 FFT 内部 remap：nqty=2
```

FFT remap 不是普通面邻居通信，可能向多个远端 rank 发送。

## 10.4 FFT 计算

一个三维 FFT 的理论复杂度：

```math
C_{\mathrm{FFT3D}}
\approx
F_{\mathrm{FFT}}
\frac{M}{R}
(\log_2N_x+\log_2N_y+\log_2N_z)
```

IK differentiation 通常包含一次正向 FFT 和三个逆向 FFT，因此：

```math
C_{\mathrm{FFT,total}}
\approx4C_{\mathrm{FFT3D}}
```

## 10.5 Poisson 和场插值

Poisson 网格计算：

```math
C_{\mathrm{Poisson}}(r)
\approx F_{\mathrm{Poisson}}M/R
```

将三个场分量插值回原子：

```math
C_{\mathrm{field}}(r)
\approx
3F_{\mathrm{field}}
N_{\mathrm{local}}(r)p^3
```

## 10.6 Grid forward

三个电场分量发送到 ghost grid：

```math
T_{\mathrm{grid-forward}}(r)
=N_{\mathrm{forward-grid}}(r)
\times3\times S_{\mathrm{FFT}}
```

double precision 时：

```math
T_{\mathrm{grid-forward}}
=24N_{\mathrm{forward-grid}}
```

## 10.7 Allreduce

PPPM 可能归约：

```text
energy：1 double
virial：6 doubles
```

当前 CCDG 的逻辑 collective 字节口径：

```math
T_{\mathrm{collective}}
=R\times count\times datatype\_bytes
```

实际 NoC link bytes 还取决于使用 ring、tree 或 recursive doubling，应在
WSE/NoC 模型中另行展开。

## 11. Reverse communication

Newton pair 开启时，ghost atom 上累积的力需要返回 owner。

当前普通原子力包含：

```text
fx、fy、fz：3 doubles/atom
```

每 rank：

```math
T_{\mathrm{reverse}}(r)
=\sum_iN_{\mathrm{reverse}}(r,i)\times3\times8
```

周期均匀体系中，平均 reverse 和 forward 数据量接近，但单个 rank 可能不同。

## 12. Output

当 timestep 命中 thermo/output 条件时，还会执行：

```text
温度、能量、压力等本地统计
全局 Reduce/Allreduce
rank 0 输出
```

输出计算量通常远小于 pair 和 PPPM，但 collective 延迟可能影响短算例。

理论模型应根据：

```text
timestep % thermo_interval
是否为 run 最后一步
```

判断该步是否包含 output。

## 13. 通用 T1 公式

稳态 short：

```math
T_1(r)
=T_{\mathrm{forward}}(r)
+T_{\mathrm{reverse}}(r)
+T_{\mathrm{pair-comm}}(r)
+T_{\mathrm{output-collective}}(r)
```

重建 short：

```math
T_1^{\mathrm{rebuild}}(r)
=T_{\mathrm{exchange}}(r)
+T_{\mathrm{borders}}(r)
+T_{\mathrm{reverse}}(r)
+T_{\mathrm{pair-comm}}(r)
+T_{\mathrm{output-collective}}(r)
```

稳态 long：

```math
T_1^{\mathrm{long}}(r)
=T_{\mathrm{CommBrick}}(r)
+T_{\mathrm{grid-reverse}}(r)
+T_{\mathrm{FFT-remap}}(r)
+T_{\mathrm{grid-forward}}(r)
+T_{\mathrm{Kspace-reduce}}(r)
+T_{\mathrm{output-collective}}(r)
```

如果目标是 NoC 链路负载：

```math
T_{\mathrm{hop}}(r)
=\sum_{m\in messages(r)}
bytes(m)\times hops(src_m,dst_m)
```

## 14. 通用 C1 公式

short：

```math
C_1(r)
=C_{\mathrm{integrate}}(r)
+C_{\mathrm{neighbor}}(r)
+C_{\mathrm{pair}}(r)
+C_{\mathrm{output}}(r)
```

EAM：

```math
C_1^{\mathrm{EAM}}(r)
=C_{\mathrm{integrate}}
+C_{\mathrm{neighbor}}
+C_{\mathrm{density}}
+C_{\mathrm{embedding}}
+C_{\mathrm{force}}
+C_{\mathrm{output}}
```

long：

```math
C_1^{\mathrm{long}}(r)
=C_{\mathrm{integrate}}
+C_{\mathrm{neighbor}}
+C_{\mathrm{real-space}}
+C_{\mathrm{rho}}
+C_{\mathrm{FFT}}
+C_{\mathrm{Poisson}}
+C_{\mathrm{field}}
+C_{\mathrm{output}}
```

## 15. LiAlOCl short、2688 原子、16 rank 示例

输入：

```text
原始 data：84 atoms
replicate：2 × 4 × 4
总原子数：2688
盒子：27.987985 × 51.971974 × 51.932017 Å
处理器网格：4 × 4 × 1
force cutoff：10.0 Å
skin：0.3 Å
ghost width：10.3 Å
```

密度：

```math
\rho
=2688/(27.987985\times51.971974\times51.932017)
\approx0.035584
```

rank 子域：

```text
lx = 6.997 Å
ly = 12.993 Å
lz = 51.932 Å
Nlocal_avg = 168
```

### 15.1 理论通信量

X 首轮：

```math
N_{x,0}=168
```

X 第二轮：

```math
N_{x,1}
\approx
\rho(10.3-6.997)l_yl_z
\approx79.31
```

Y 首轮：

```math
N_{y,0}
\approx
\rho(10.3)(l_x+2\times10.3)l_z
\approx525.28
```

每 rank 发送 atom-copies：

```math
N_{\mathrm{copies}}
\approx2(168+79.31+525.28)
\approx1545.16
```

Forward：

```math
T_{\mathrm{forward}}
\approx1545.16\times3\times8
\approx37084\ bytes
```

Reverse：

```math
T_{\mathrm{reverse}}
\approx37084\ bytes
```

总 T1：

```math
T_1\approx74168\ bytes/rank/timestep
```

实际 WSE plan：

```text
全局发送量：1,175,040 bytes
平均每 rank：73,440 bytes
理论误差：约 0.99%
```

### 15.2 理论计算量

Neighbor list cutoff 为 10.3 Å：

```math
N_{\mathrm{list}}
\approx
\frac12\times168\times
0.035584\times\frac{4\pi}{3}(10.3)^3
\approx13681/rank
```

实际日志：

```text
平均 13,544/rank
最少 13,044/rank
最多 14,286/rank
```

Force cutoff 为 10 Å：

```math
N_{\mathrm{pair}}
\approx
\frac12\times168\times
0.035584\times\frac{4\pi}{3}(10.0)^3
\approx12520/rank
```

Pair 模型：

```math
C_{\mathrm{pair}}
\approx8\times13681
+(45\sim50)\times12520
```

得到：

```text
Cpair ≈ 0.673～0.735 M scalar operations/rank
```

积分：

```math
C_{\mathrm{integrate}}
\approx18\times168
=3024
```

该 timestep 没有 neighbor rebuild，因此：

```text
Cneighbor = 0
```

最终：

```text
C1 ≈ 0.676～0.738 M scalar operations/rank/timestep
```

更稳健的工作量表达：

```text
Nlocal ≈ 168/rank
neighbor entries ≈ 13,681/rank
force pair evaluations ≈ 12,520/rank
```

## 16. 理论值与运行值如何对齐

通信量 T1：

```text
理论 estimator
  → 按 rank、phase、dst 输出 messages/bytes

wse_plan.json
  → 读取 scope=run 的真实 pack 字节数

对比：
  每 rank
  每 phase
  每 destination
  总 bytes
```

建议验收：

```text
每 rank、每 phase、每 dst 的字节误差 ≤ 2%
```

计算量 C1：

```text
理论 estimator
  → pair evaluations、grid updates、FFT work、modeled ops

运行侧
  → LAMMPS Pair/Neigh/Kspace timer
  → PAPI/LIKWID/perf 浮点和 cycles 计数
```

建议验收：

```text
neighbor/pair 数误差 ≤ 5%
校准后计算周期或时间误差 ≤ 10%
```

不能直接用 DUMPI 的 MPI 调用间隔作为真实 FLOP 数；它更适合验证时间线和通信。

## 17. 估算器需要输出的最小数据

```json
{
  "rank": 0,
  "nlocal": 168,
  "nghost_est": 2184,
  "communication": {
    "forward_bytes": 37084,
    "reverse_bytes": 37084,
    "pair_comm_bytes": 0,
    "grid_reverse_bytes": 0,
    "fft_remap_bytes": 0,
    "grid_forward_bytes": 0,
    "collective_bytes": 0,
    "total_send_bytes": 74168
  },
  "computation": {
    "neighbor_entries": 13681,
    "pair_evaluations": 12520,
    "integrate_ops": 3024,
    "modeled_ops_min": 675872,
    "modeled_ops_max": 738472
  }
}
```

同时必须保存估算假设：

```json
{
  "scope": "run",
  "step_type": "steady",
  "uniform_density": true,
  "newton_pair": true,
  "send_only": true,
  "force_cutoff": 10.0,
  "neighbor_cutoff": 10.3,
  "sqrt_div_cost_model": "one scalar operation each"
}
```

只有明确这些假设，T1/C1 才能在不同实验和不同实现之间进行比较。
