# attic/ — 历史/一次性实验脚本（未纳入复现路径）

这些脚本硬编码了原始工作区的历史 run 目录（未随包分发），仅作方法存档：
- run_wse_main.sh / run_wse_traffic.sh / cap_scan_free.sh: 早期 WSE 门控实验
- run_booksim_batch.sh: 旧批量注入（已被 run_noc_pipeline.sh 取代）
- run_kspace_profile.sh: SMPI mesh 实验（需 lammps-build-smpi，未随包）
- gen_marching_ccdg.py: marching CCDG 生成器原型
- validate_simgrid.py: 旧批量 SimGrid 验证（已被 ccdg_simgrid_validate.py 取代）

主流水线请用 ../run_noc_pipeline.sh。
