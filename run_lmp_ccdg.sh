#!/usr/bin/env bash
# ============================================================================
# LAMMPS → DUMPI trace → CCDG → SimGrid 验证 流水线引擎
# （run_short_lmp.sh / run_long_lmp.sh 的公共实现，请勿直接调用）
#
# 直接用法（推荐走包装脚本）:
#   ./run_short_lmp.sh <rank数> <cu|h2o|lialocl> <原子数> [BookSim步数=10] [sim=none]   # 无长程力 Kspace
#   ./run_long_lmp.sh  <rank数> <cu|h2o|lialocl> <原子数> [BookSim步数=10] [sim=none]   # 含长程力 Kspace(PPPM)
#   sim=none|booksim|demand: demand = 需求驱动数据流编译（Cerebras 范式）+
#   BookSim 两档对比（trace free 基线 vs 无同步 demand 档）
#
# 流程（真实 LAMMPS 只跑一次、固定 1 步，以控制 trace/CCDG 的 MPI 节点规模）:
#   ① 按 体系/模式/原子数 生成 in.lammps（捕获用 run 1，可用 CAPTURE_STEPS 环境变量覆盖）
#   ② 原生 mpirun(真实 MPI) 运行 LAMMPS 这一次，DUMPI(LD_PRELOAD) 截获通信 trace
#   ③ dumpi2ccdg 生成原始 CCDG + 压缩 CCDG（setup 裁剪 + 突发折叠, 仅含迭代段）
#      + trimonly CCDG（setup 裁剪、不折叠，需求驱动编译器的输入载体）
#   ④ SimGrid DAG 仿真验证压缩版单步 CCDG: |T_sim - T_real| / T_real ≤ 5% 记为 PASS
#      (T_real 取 Loop time，纯迭代段，与裁剪后 CCDG 口径对齐)
#   ⑤ [sim=booksim] 压缩单步 CCDG 用 ccdg_unroll.py 展开到 BookSim步数，注入 BookSim
#      方形 2D mesh 测 total_cycles 与平均每轮迭代 cycles
#      (rank 须为完全平方数 4/9/16/25; 超时秒数用 BOOKSIM_TIMEOUT 环境变量, 默认 5400)
#   ⑥ [sim=demand] trimonly → ccdg_demand.py 需求提炼 + 时槽排布 → 无同步 CCDG(+est 表)，
#      BookSim 两档对比: trace free 基线 vs demand 档（blocked=0 by construction）
#
# 退出码: 0=PASS  1=FAIL(误差>5%)  2=流水线错误(缺文件/运行失败)
# 产物:   runs/pipeline/<mode>_<体系>_<原子数>a_<rank数>r_<时间戳>/
#         ├── in.lammps / lammps.log / dumpi-*.bin|meta   (单步捕获)
#         ├── trace_<rank数>ranks_global.ccdg              (原始，含 setup)
#         ├── compact_<rank数>ranks_global.ccdg            (裁剪+压缩，仅迭代段)
#         ├── trimonly_<rank数>ranks_global.ccdg           (裁剪不折叠，仅迭代段)
#         ├── simgrid_dag/{platform.xml, ccdg_dag_sim.cpp, ccdg_dag_sim}
#         ├── validation_result.txt
#         ├── unrolled_<BookSim步数>steps.ccdg             # [sim=booksim]
#         ├── booksim_result.txt                            # [sim=booksim]
#         ├── demand_plan.log                               # [sim=demand] 需求提炼+排布报告
#         ├── booksim_free/{trimonly_*.ccdg}                # [sim=demand] free 基线档
#         └── booksim_demand/{demand_*.ccdg, demand_*.est}  # [sim=demand] 无同步 demand 档
#
# 说明:
#   - cu 用 EAM 势，本身无长程库仑项(Kspace 不适用)，short/long 两模式结果一致
#   - 为控制 DUMPI trace 规模，生成输入不含 minimize，仅 run 段(NVE)
# ============================================================================
set -uo pipefail

ROOT="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
INSTALL=$ROOT/install
LMP=$INSTALL/bin/lmp
LIBDUMPI=$INSTALL/lib/libdumpi.so
DUMPI2CCDG=$ROOT/dumpi2ccdg/dumpi2ccdg
VALIDATOR=$ROOT/ccdg_simgrid_validate.py
RUNS_ROOT=$ROOT/runs/pipeline

MODE=${1:-}
RANKS=${2:-}
SYSTEM=${3:-}
NATOMS=${4:-}
STEPS=${5:-10}          # BookSim 仿真步数（仅 sim=booksim 时用于展开单步压缩 CCDG）
SIM=${6:-none}
CAPTURE_STEPS=${CAPTURE_STEPS:-1}   # 真实 LAMMPS 捕获/验证步数，固定跑这一次

usage() {
  sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
}

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# ── 参数校验 ──────────────────────────────────────────────────────────────
case "$MODE" in short|long) ;; *) usage ;; esac
[[ "$RANKS" =~ ^[0-9]+$ ]] && [ "$RANKS" -ge 1 ] || { echo "ERROR: rank数须为正整数"; exit 2; }
case "$SYSTEM" in cu|h2o|lialocl) ;; *) echo "ERROR: 体系须为 cu|h2o|lialocl"; usage ;; esac
[[ "$NATOMS" =~ ^[0-9]+$ ]] && [ "$NATOMS" -ge 1 ] || { echo "ERROR: 原子数须为正整数"; exit 2; }
[[ "$STEPS" =~ ^[0-9]+$ ]] && [ "$STEPS" -ge 1 ] || { echo "ERROR: BookSim 步数须为正整数"; exit 2; }
[[ "$CAPTURE_STEPS" =~ ^[0-9]+$ ]] && [ "$CAPTURE_STEPS" -ge 1 ] || { echo "ERROR: CAPTURE_STEPS 须为正整数"; exit 2; }
case "$SIM" in none|booksim|demand) ;; *) echo "ERROR: sim 参数须为 none|booksim|demand"; usage ;; esac

for f in "$LMP" "$LIBDUMPI" "$DUMPI2CCDG" "$VALIDATOR"; do
  [ -e "$f" ] || { echo "ERROR: 缺少组件 $f"; exit 2; }
done

PAIR_SUFFIX=$MODE
[ "$MODE" = short ] && PAIR_SUFFIX=cut
[ "$MODE" = long ]  && PAIR_SUFFIX=long

RUN_DIR=$RUNS_ROOT/${MODE}_${SYSTEM}_${NATOMS}a_${RANKS}r_$(date +%Y%m%d_%H%M%S)
mkdir -p "$RUN_DIR"
log "模式=$MODE 体系=$SYSTEM rank=$RANKS 目标原子=$NATOMS 捕获步数=$CAPTURE_STEPS BookSim步数=$STEPS sim=$SIM"
log "运行目录: $RUN_DIR"

# ── ① 生成 in.lammps ─────────────────────────────────────────────────────
# make_input <步数> <输出目录>
make_input() {
  local steps=$1 outdir=$2
  case "$SYSTEM" in
  cu)
    # fcc N^3×4 ≈ NATOMS
    local n actual
    n=$(awk -v t="$NATOMS" 'BEGIN{n=int((t/4)^(1/3)+0.5); if(n<1)n=1; printf "%d",n}')
    actual=$((4*n*n*n))
    [ "$MODE" = long ] && log "注意: Cu/EAM 无长程力(Kspace 不适用)，long 与 short 结果一致"
    cp "$ROOT/cases/eam_cu_1k/Cu_u3.eam" "$outdir/"
    cat > "$outdir/in.lammps" <<EOF
# Cu EAM benchmark (short/long 相同): target=${NATOMS}, actual=${actual} (fcc ${n}^3*4)
units           metal
boundary        p p p
atom_style      atomic

lattice         fcc 3.615
region          box block 0 ${n} 0 ${n} 0 ${n}
create_box      1 box
create_atoms    1 box

mass            1 63.546
pair_style      eam
pair_coeff      * * Cu_u3.eam

velocity        all create 300.0 12345 units box

timestep        0.001
neighbor        0.3 bin
neigh_modify    delay 0 every 1
fix             1 all nve

thermo          10
run             ${steps}
EOF
    ;;
  h2o)
    local n actual cut
    n=$(awk -v t="$NATOMS" 'BEGIN{n=int((t/4)^(1/3)+0.5); if(n<1)n=1; printf "%d",n}')
    actual=$((4*n*n*n))
    if [ $((actual % 2)) -eq 1 ]; then actual=$((actual-1)); fi  # 保证净电荷为 0(偶原子数)
    cut=$(awk -v n="$n" 'BEGIN{c=2*n-0.3; if(c>10)c=10; printf "%.2f",c}')
    cat > "$outdir/in.lammps" <<EOF
# H2O 类 2-型 LJ+Coul 体系: target=${NATOMS}, actual=${actual} (fcc ${n}^3*4), mode=${MODE}
units           metal
boundary        p p p
atom_style      charge

lattice         fcc 4.0
region          box block 0 ${n} 0 ${n} 0 ${n}
create_box      2 box
create_atoms    1 box
group           type2 id 2:${actual}:2
set             group type2 type 2

set             type 1 charge 0.01
set             type 2 charge -0.01
mass            1 1.0
mass            2 1.0

pair_style      lj/cut/coul/${PAIR_SUFFIX} ${cut}${MODE:0:0}$([ "$MODE" = long ] && echo " ${cut}")
pair_coeff      1 1 1.0 1.0 ${cut}
pair_coeff      2 2 1.0 1.0 ${cut}
pair_coeff      1 2 1.0 1.0 ${cut}
$([ "$MODE" = long ] && echo "kspace_style    pppm 1.0e-4")

timestep        0.001
velocity        all create 0.1 12345

neighbor        0.3 bin
neigh_modify    delay 0 every 1
fix             1 all nve

thermo          10
run             ${steps}
EOF
    ;;
  lialocl)
    local a b c actual
    read -r a b c actual <<<"$(awk -v t="$NATOMS" 'BEGIN{
      n=int((t/84)^(1/3))+1; if(n<1)n=1; best=1e18;
      for(x=1;x<=n+2;x++) for(y=1;y<=n+2;y++) for(z=1;z<=n+2;z++){
        d=84*x*y*z; df=d-t; if(df<0) df=-df; if(df<best){best=df;ba=x;bb=y;bc=z}}
      printf "%d %d %d %d", ba, bb, bc, 84*ba*bb*bc}')"
    cp "$ROOT/cases/lialocl_coul/data.LiAlOCl_nvt_charge" "$outdir/"
    cat > "$outdir/in.lammps" <<EOF
# LiAlOCl LJ+Coul 4-型体系: target=${NATOMS}, actual=${actual} (replicate ${a} ${b} ${c}), mode=${MODE}
units           metal
boundary        p p p
atom_style      charge

read_data       ./data.LiAlOCl_nvt_charge
replicate       ${a} ${b} ${c}

set             type 1 charge 1.0
set             type 2 charge 3.0
set             type 3 charge -1.0
set             type 4 charge -2.0

pair_style      lj/cut/coul/${PAIR_SUFFIX} 10.0$([ "$MODE" = long ] && echo " 10.0")
pair_coeff      1 1 0.001 1.506
pair_coeff      2 2 0.010 1.044
pair_coeff      3 3 0.100 4.417
pair_coeff      4 4 0.150 3.500
$([ "$MODE" = long ] && echo "kspace_style    pppm 1.0e-4")

timestep        0.001
velocity        all create 600.0 12345

neighbor        0.3 bin
neigh_modify    delay 0 every 1
fix             1 all nve

thermo          10
run             ${steps}
EOF
    ;;
  esac
  echo "${actual}"
}
ACTUAL_ATOMS=$(make_input "$CAPTURE_STEPS" "$RUN_DIR")
log "in.lammps 已生成 (实际原子 ${ACTUAL_ATOMS})"

# ── ② DUMPI capture + 真实 LAMMPS 运行 ───────────────────────────────────
export LD_PRELOAD=$LIBDUMPI
export LD_LIBRARY_PATH=$INSTALL/lib:${LD_LIBRARY_PATH:-}
export DUMPI_OUTDIR=$RUN_DIR

log "运行真实 LAMMPS (mpirun -np $RANKS, DUMPI 截获中)..."
( cd "$RUN_DIR" && timeout 1800 mpirun -np "$RANKS" --allow-run-as-root --oversubscribe \
    -x LD_PRELOAD -x LD_LIBRARY_PATH -x DUMPI_OUTDIR \
    "$LMP" -in in.lammps > lammps.log 2>&1 )
RC=$?
if [ $RC -ne 0 ]; then
  echo "ERROR: LAMMPS 退出码 $RC，日志尾部:" >&2
  tail -20 "$RUN_DIR/lammps.log" >&2
  exit 2
fi
NMETA=$(ls "$RUN_DIR"/dumpi-*.meta 2>/dev/null | wc -l)
echo "  dumpi-*.meta 文件数: $NMETA"
[ "$NMETA" -ge 1 ] || { echo "ERROR: 未捕获到 DUMPI trace"; exit 2; }
grep -q "Loop time" "$RUN_DIR/lammps.log" || { echo "ERROR: log 中无 Loop time"; exit 2; }
grep "Loop time" "$RUN_DIR/lammps.log"

# ── ③ dumpi2ccdg ─────────────────────────────────────────────────────────
CCDG=$RUN_DIR/trace_${RANKS}ranks_global.ccdg
log "生成 CCDG ..."
LD_LIBRARY_PATH=$INSTALL/lib:${LD_LIBRARY_PATH:-} \
  "$DUMPI2CCDG" "$RUN_DIR" > "$CCDG" 2> "$RUN_DIR/ccdg_gen.log"
grep -q '"num_ranks"' "$CCDG" || {
  echo "ERROR: CCDG 生成失败"; tail -5 "$RUN_DIR/ccdg_gen.log"; exit 2; }
echo "  CCDG: $CCDG ($(wc -c < "$CCDG") bytes, ranks=$(grep -o '"num_ranks":[0-9]*' "$CCDG" | head -1))"

# ── ③b 压缩 CCDG：setup 裁剪 + 突发折叠（仅迭代段，供 SimGrid/BookSim 使用）──
CCDG_COMPACT=$RUN_DIR/compact_${RANKS}ranks_global.ccdg
log "生成压缩 CCDG (CCDG_TRIM_SETUP + CCDG_COMPACT) ..."
# 注意 dumpi2ccdg 的 getenv 只查存在性不查值：CCDG_COMPACT 赋值任意值（含 0）都触发折叠
CCDG_TRIM_SETUP=1 CCDG_COMPACT=1 LD_LIBRARY_PATH=$INSTALL/lib:${LD_LIBRARY_PATH:-} \
  "$DUMPI2CCDG" "$RUN_DIR" > "$CCDG_COMPACT" 2> "$RUN_DIR/ccdg_compact_gen.log"
grep -q '"num_ranks"' "$CCDG_COMPACT" || {
  echo "ERROR: 压缩 CCDG 生成失败"; tail -5 "$RUN_DIR/ccdg_compact_gen.log"; exit 2; }
grep -E "Trim setup: removed|Compact bursts" "$RUN_DIR/ccdg_compact_gen.log" | sed 's/^/  /'
echo "  压缩 CCDG: $CCDG_COMPACT ($(wc -c < "$CCDG_COMPACT") bytes)"

# ── ③c trimonly CCDG：setup 裁剪、不折叠（需求驱动编译器的输入载体）──
# 必须 env -u CCDG_COMPACT 清除外部残留变量：getenv 只查存在性，置 0 仍会折叠
CCDG_TRIMONLY=$RUN_DIR/trimonly_${RANKS}ranks_global.ccdg
log "生成 trimonly CCDG (CCDG_TRIM_SETUP, 无 COMPACT) ..."
env -u CCDG_COMPACT CCDG_TRIM_SETUP=1 LD_LIBRARY_PATH=$INSTALL/lib:${LD_LIBRARY_PATH:-} \
  "$DUMPI2CCDG" "$RUN_DIR" > "$CCDG_TRIMONLY" 2> "$RUN_DIR/ccdg_trimonly_gen.log"
grep -q '"num_ranks"' "$CCDG_TRIMONLY" || {
  echo "ERROR: trimonly CCDG 生成失败"; tail -5 "$RUN_DIR/ccdg_trimonly_gen.log"; exit 2; }
echo "  trimonly CCDG: $CCDG_TRIMONLY ($(wc -c < "$CCDG_TRIMONLY") bytes)"

# ── ④ SimGrid DAG 验证（压缩版单步 CCDG vs Loop time）────────────────────────────
# 未设置 SIMGRID_CPU_FREQ 时自动迭代校准频率；设置则用固定频率模式
if [ -n "${SIMGRID_CPU_FREQ:-}" ]; then
  export SIMGRID_CPU_FREQ
  log "SimGrid DAG 验证 (固定 CPU=${SIMGRID_CPU_FREQ} GHz) ..."
else
  log "SimGrid DAG 验证 (自动校准频率) ..."
fi
python3 "$VALIDATOR" "$RUN_DIR" "$RANKS" "$CCDG_COMPACT" | tee "$RUN_DIR/validation_result.txt"
RC=${PIPESTATUS[0]}

echo ""
echo "==== 完成 ===="
echo "结果: $RUN_DIR/validation_result.txt"
grep -E "T_real|T_simgrid|E_total" "$RUN_DIR/validation_result.txt" || true
case $RC in
  0) log "FINAL: PASS (误差 ≤ 5%)"; ;;
  1) log "FINAL: FAIL (误差 > 5%)"; ;;
  *) log "FINAL: ERROR"; exit 2; ;;
esac

# ── ⑤ [可选] BookSim 方形 2D mesh NoC 仿真 ─────────────────────────────
# 口径: cpu=noc=2.0GHz, flit_size_bytes=1; 结果追加 booksim2/results/ccdg_mesh_results.csv
# 压缩单步 CCDG 用 ccdg_unroll.py 展开到 STEPS 步（无需重跑 LAMMPS），
# 注入后统计 total_cycles 与平均每轮迭代 cycles (cycles_per_iter)
if [ "$SIM" = booksim ]; then
  BS_RUNNER=$ROOT/booksim2/run_ccdg_mesh.sh
  UNROLLER=$ROOT/ccdg_unroll.py
  BS_TIMEOUT=${BOOKSIM_TIMEOUT:-5400}
  if [ ! -x "$BS_RUNNER" ] || [ ! -f "$UNROLLER" ]; then
    echo "ERROR: 缺少 BookSim 注入脚本 $BS_RUNNER 或展开器 $UNROLLER" >&2
    exit 2
  fi

  UNROLLED=$RUN_DIR/unrolled_${STEPS}steps.ccdg
  log "展开压缩单步 CCDG 到 $STEPS 步 ..."
  python3 "$UNROLLER" "$CCDG_COMPACT" "$STEPS" "$UNROLLED" || exit 2

  log "BookSim mesh 仿真 (超时 ${BS_TIMEOUT}s, 大规模 CCDG 可能需小时级) ..."
  BS_RC=0
  "$BS_RUNNER" "$UNROLLED" "$BS_TIMEOUT" "$STEPS" 2>&1 | tee "$RUN_DIR/booksim_result.txt" || BS_RC=${PIPESTATUS[0]}
  BS_CSV=$ROOT/booksim2/results/ccdg_mesh_results.csv
  if tail -1 "$BS_CSV" 2>/dev/null | grep -q "$(basename "$RUN_DIR")"; then
    log "BookSim: $(tail -1 "$BS_CSV")"
  else
    log "BookSim: 仿真未产出结果 (rc=$BS_RC)，详见 $RUN_DIR/booksim_result.txt"
  fi
fi

# ── ⑥ [可选] 需求驱动数据流编译（Cerebras 范式）+ BookSim 两档对比 ──────────
# 口径: trimonly 单步 trace free 基线 vs 无同步 demand 档（各注入 1 步，STEPS 不适用）。
# 两档 CCDG 放不同子目录: run_ccdg_mesh.sh 的 TAG 取父目录名，同目录会互相覆盖。
if [ "$SIM" = demand ]; then
  DEMAND_CC=$ROOT/booksim2/ccdg_demand.py
  BS_RUNNER=$ROOT/booksim2/run_ccdg_mesh.sh
  BS_TIMEOUT=${BOOKSIM_TIMEOUT:-5400}
  if [ ! -f "$DEMAND_CC" ] || [ ! -x "$BS_RUNNER" ]; then
    echo "ERROR: 缺少需求驱动编译器 $DEMAND_CC 或注入脚本 $BS_RUNNER" >&2
    exit 2
  fi

  FREE_DIR=$RUN_DIR/booksim_free
  DEMAND_DIR=$RUN_DIR/booksim_demand
  mkdir -p "$FREE_DIR" "$DEMAND_DIR"
  cp "$CCDG_TRIMONLY" "$FREE_DIR/trimonly_${RANKS}ranks_global.ccdg"

  log "需求提炼 + 时槽排布 + 无同步 CCDG 生成 (ccdg_demand.py) ..."
  python3 "$DEMAND_CC" "$CCDG_TRIMONLY" \
    -o "$DEMAND_DIR/demand_${RANKS}ranks" \
    --cap "${CCDG_COMPUTE_CAP:-2.5e10}" 2>&1 | tee "$RUN_DIR/demand_plan.log" || exit 2
  DEMAND_CCDG=$DEMAND_DIR/demand_${RANKS}ranks_demand.ccdg
  DEMAND_EST=$DEMAND_DIR/demand_${RANKS}ranks_demand.est
  grep -q "plan report" "$RUN_DIR/demand_plan.log" && \
    [ -f "$DEMAND_CCDG" ] && [ -f "$DEMAND_EST" ] || {
    echo "ERROR: demand CCDG 产物缺失"; exit 2; }

  log "BookSim: trace free 基线 (trimonly 1 步) ..."
  "$BS_RUNNER" "$FREE_DIR/trimonly_${RANKS}ranks_global.ccdg" "$BS_TIMEOUT" 1 \
    2>&1 | tee "$RUN_DIR/booksim_free_result.txt" || true

  log "BookSim: demand 档 (无同步 CCDG + est 门控, 1 步) ..."
  CCDG_SCHED_FILE="$DEMAND_EST" "$BS_RUNNER" "$DEMAND_CCDG" "$BS_TIMEOUT" 1 \
    2>&1 | tee "$RUN_DIR/booksim_demand_result.txt" || true

  log "BookSim 结果 (CSV 末 2 行):"
  tail -2 "$ROOT/booksim2/results/ccdg_mesh_results.csv"
fi
exit $RC