#!/usr/bin/env bash
# 含长程力 Kspace(PPPM) 的 LAMMPS → DUMPI → CCDG → SimGrid 验证流水线（可选 BookSim NoC 仿真）
# 用法: ./run_long_lmp.sh <rank数> <cu|h2o|lialocl> <原子数> [BookSim步数=10] [sim=none|booksim|demand]
# 说明: 真实 LAMMPS 固定只跑 1 步，SimGrid 验证裁剪 setup+压缩后的单步 CCDG（纯迭代段）;
#       步数参数仅 sim=booksim 时生效——压缩单步 CCDG 展开为该步数注入 BookSim，
#       统计 total_cycles 与平均每轮迭代 cycles（无需重跑 LAMMPS）;
#       sim=demand 时额外生成 trimonly CCDG 并走需求驱动数据流编译（ccdg_demand.py），
#       BookSim 注入两档对比: trace free 基线 vs 无同步 demand 档（各 1 步）
# 例:   ./run_long_lmp.sh 64 lialocl 2688                 # 单步捕获 + SimGrid 验证
#       ./run_long_lmp.sh 25 lialocl 2688 100 booksim     # 验证单步 + BookSim 展开 100 步仿真
#       ./run_long_lmp.sh 25 lialocl 2688 1 demand        # + 需求驱动编译两档 BookSim 对比
# 注意: cu 用 EAM 势，无长程库仑项，与 short 模式结果一致
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_lmp_ccdg.sh" long "$@"