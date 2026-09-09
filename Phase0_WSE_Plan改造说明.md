# Phase 0：LAMMPS WSE Plan 改造说明

## 1. 改造目标

原流水线只能通过 DUMPI 截获 MPI 调用，再转换成 CCDG：

```text
in.lammps → LAMMPS + DUMPI → dumpi-*.bin/meta
          → dumpi2ccdg → trimonly.ccdg → demand 编译 / BookSim
```

Phase 0 在 LAMMPS 知道通信对象和数据量的位置直接导出通信计划：

```text
in.lammps → LAMMPS
              ├→ CommBrick rank 分片
              ├→ Kspace/FFT rank 分片
              └→ merge_wse_plan.py → wse_plan.json

           LAMMPS + DUMPI → trimonly.ccdg
                                  │
                    validate_wse_plan.py
                                  │
                    对拍 wse_plan.json
```

当前 DUMPI/CCDG 仍是 BookSim 的输入，也是 `wse_plan.json` 的校准 oracle。
Phase 0 尚未改变后续编译器输入；Phase 1 才会让 WSE 编译器直接读取
`wse_plan.json`。

## 2. 设计原则

1. **不增加 MPI 通信**：只读取 LAMMPS 已经计算出的目标 rank、count 和
   datatype 大小。
2. **每个 rank 独立写文件**：避免共享文件锁和 rank 间同步。
3. **区分 setup 和 run**：初始化通信保留在 plan 中，但对拍和编译主要使用
   `scope="run"`。
4. **以实际 pack 结果计算字节数**：不根据原子数猜测消息大小。
5. **DUMPI 降级为 oracle**：源码 plan 与 trim-only CCDG 不一致时立即报错。
6. **只对齐现有脚本场景**：支持 `K×K×1`、CommBrick、Cu/EAM，以及
   H2O/LiAlOCl 的 short 和 PPPM long 路径。

## 3. LAMMPS 源码改造

### 3.1 CommBrick 通信

涉及文件：

- `lammps-src/src/comm_brick.h`
- `lammps-src/src/comm_brick.cpp`

通过环境变量开启：

```bash
export LAMMPS_WSE_PLAN=/path/to/run/wse_plan
```

每个 rank 生成：

```text
wse_plan.rank0000.jsonl
wse_plan.rank0001.jsonl
...
```

当前记录：

- `swap_setup`：处理器网格、swap 维度、方向、轮次和 PBC 信息；
- `borders`：建立 ghost atom 时的消息；
- `forward`：坐标/属性正向通信；
- `reverse`：力的反向归并；
- `pair_forward`、`pair_reverse`：Cu/EAM 等 Pair 样式的附加通信。

### 3.2 PPPM Grid3d 与 FFT

涉及文件：

- `lammps-src/src/wse_plan_kspace.h`
- `lammps-src/src/grid3d.cpp`
- `lammps-src/src/KSPACE/remap.cpp`
- `lammps-src/src/KSPACE/pppm.cpp`

Kspace 使用独立分片：

```text
wse_plan.kspace.rank0000.jsonl
wse_plan.kspace.rank0001.jsonl
...
```

使用独立分片是为了避免 CommBrick 和 Kspace 的两个带缓冲文件句柄同时写同一
文件。离线合并后仍只有一个 `wse_plan.json`。

当前记录：

- `grid_reverse`：PPPM 电荷密度 ghost grid 向 owner 汇总；
- `grid_forward`：电场网格向 ghost 区域传播；
- `fft_remap`：`remap_3d()` 中全部 FFT transpose 点对点消息；
- `kspace_reduce`：PPPM energy/virial Allreduce。

`remap_3d()` 是统一插桩点，因为 `brick2fft` 和 FFT 内部多次 transpose 最终都会
经过这里。只修改上层 `Remap::perform()` 会遗漏 FFT 内部通信。

### 3.3 Kspace writer

`WsePlanKspace` 是 header-only singleton，负责：

- 保存当前 rank、timestep 和 setup/run 上下文；
- 让 `pppm.cpp`、`grid3d.cpp` 和 `remap.cpp` 共用同一文件句柄；
- 输出 message 和 collective 记录；
- 在进程结束时关闭文件。

它不调用 `MPI_Comm_rank` 或额外的 collective，rank 信息由 `PPPM::compute()`
直接传入。

## 4. Plan 数据格式

### 4.1 Metadata

```json
{
  "kind": "metadata",
  "component": "kspace",
  "schema_version": 1,
  "rank": 0,
  "num_ranks": 16
}
```

CommBrick metadata 还包含 `procgrid` 和 `myloc`。

### 4.2 点对点消息

```json
{
  "kind": "message",
  "component": "kspace",
  "scope": "run",
  "phase": "fft_remap",
  "timestep": 1,
  "src": 0,
  "dst": 5,
  "value_count": 162,
  "datatype_bytes": 8,
  "bytes": 1296
}
```

计算关系：

```text
bytes = value_count × datatype_bytes
```

### 4.3 Collective

```json
{
  "kind": "collective",
  "component": "kspace",
  "scope": "run",
  "phase": "kspace_reduce",
  "operation": "allreduce",
  "value_count": 6,
  "datatype_bytes": 8,
  "bytes": 768
}
```

为与当前 CCDG 口径一致：

```text
bytes = num_ranks × value_count × datatype_bytes
```

## 5. 脚本链路变化

`run_noc_pipeline.sh` 默认设置：

```bash
WSE_PLAN_CAPTURE=1
```

如需关闭：

```bash
WSE_PLAN_CAPTURE=0 ./run_noc_pipeline.sh ...
```

运行阶段新增三个步骤：

1. 将 `LAMMPS_WSE_PLAN` 传递给全部 MPI rank；
2. 调用 `merge_wse_plan.py` 合并 CommBrick 和 Kspace 分片；
3. 调用 `validate_wse_plan.py` 与 `trimonly_*.ccdg` 对拍。

新增产物：

```text
run_dir/
├── wse_plan.rankNNNN.jsonl
├── wse_plan.kspace.rankNNNN.jsonl   # long/PPPM
├── wse_plan.json
├── wse_plan_merge.log
└── wse_plan_validation.json
```

`wse_plan.json` 包含：

- rank 和处理器网格 metadata；
- CommBrick 与 Kspace component metadata；
- setup/run 全部记录；
- 按 scope、phase 汇总的消息数、collective 数和字节数。

## 6. 对齐方法

`validate_wse_plan.py` 只比较 `scope="run"`：

### 点对点通信

1. 从 plan 读取所有 message；
2. 从 trim-only CCDG 读取 `SEND/ISEND`；
3. 按 `(src,dst)` 映射成二维方向；
4. 比较每个方向的消息数和字节数；
5. 每个方向字节偏差不超过 2% 才通过。

### Kspace collective

1. 从 plan 读取 Kspace collective 的操作和字节数；
2. 检查 CCDG 中存在足够数量的同类型、同字节数 collective；
3. thermo 等非 Kspace collective 可以额外存在，不会误判为 Kspace。

## 7. 已完成的对齐实验

以下命令均完成 LAMMPS、DUMPI、CCDG、WSE plan 对拍和 BookSim free 仿真：

```bash
./run_noc_pipeline.sh short 16 lialocl 2688 free
./run_noc_pipeline.sh short 4  h2o     256 free
./run_noc_pipeline.sh short 4  cu      256 free
./run_noc_pipeline.sh long  4  h2o     256 free
./run_noc_pipeline.sh long  4  lialocl 84   free
./run_noc_pipeline.sh long  16 lialocl 2688 free
```

关键结果：

- short 16-rank LiAlOCl：192 条消息，1,175,040 bytes；
- short 4-rank Cu/EAM：64 条消息，66,560 bytes；
- long 4-rank H2O：158 条消息，167,208 bytes；
- long 4-rank LiAlOCl：158 条消息，292,240 bytes；
- long 16-rank LiAlOCl：1,220 条消息，7,823,808 bytes。

上述实验的点对点消息数和逐方向字节数均与 trim-only CCDG 完全一致；
PPPM energy/virial collective 也被 CCDG 覆盖。short 回归与 BookSim unresolved
检查均通过。

## 8. 当前边界

当前实现针对现有 pipeline，不是所有 LAMMPS 通信的通用拦截层：

- DUMPI/CCDG 暂未从主流程删除；
- thermo、通用 Fix/Bond 和其他 Kspace style 尚未全部源码化；
- atom migration 的非零 exchange/Irregular 特殊路径尚未单独完善；
- `wse_plan.json` 尚未直接输入 `ccdg_demand.py` 或 BookSim；
- CommBrick 和 Kspace 暂时维护两个 writer，后续可抽取统一 writer。

因此，当前 Phase 0 的准确表述是：**现有 short 与 PPPM long 实验场景的源码级
通信计划已经形成并通过 trace 对拍，下一步可在此格式上开发 Phase 1 编译器。**
