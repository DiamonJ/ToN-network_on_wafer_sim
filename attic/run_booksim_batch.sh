cd /work1/jiangtao/lammps_trace/booksim2
CCDG=/work1/jiangtao/lammps_trace/runs/pipeline/short_lialocl_2688a_256r_20260827_153703/unrolled_1steps.ccdg

# 1) free 基线（wse_gating=0，可省略后三参数）
./run_ccdg_mesh.sh $CCDG 5400 1 0

# 2) WSE pw 扫描（wse_gating=1, phase_width=pw, strip_width=2）
for pw in 1 2 4; do
  ./run_ccdg_mesh.sh $CCDG 5400 1 1 $pw 2
done