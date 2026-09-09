#!/usr/bin/env bash
# ============================================================================
# run_noc_pipeline.sh — LAMMPS → DUMPI trace → CCDG → (free|demand) → BookSim → evaluation
#
# 无 SimGrid 版流水线：数据质量闸门改用 dumpi2ccdg 裁剪日志的 BARRIER 锚定证据
# （verlet run 首尾 barrier 对 vs log 的 Loop time）+ 注入后 unresolved/sent==recv 硬校验。
#
# 用法:
#   ./run_noc_pipeline.sh <short|long> <rank数> <cu|h2o|lialocl> <原子数> [模式=both]
#   模式: free   = trace 自由注入（trimonly 载体，含同步节点，反应式上界）
#         demand = 三档需求驱动编译: hb(正确性保证) / cerebras(XY stage+wavelet)
#                  / ilv(相位交错对照)
#         both   = 同源四档对比（默认；跨档比较只在同一次捕获内有意义）
#   载体: Z 压缩口径（processors K K 1，Z 邻居本地化），闸门校验 grid 与 b≤2
#
# 环境变量:
#   CAPTURE_STEPS     真实 LAMMPS 捕获步数（默认 1；两档注入按 1 个捕获步口径）
#   BOOKSIM_TIMEOUT   BookSim 单档超时秒（默认 5400）
#   CCDG_COMPUTE_CAP  计算能力 ops/s（默认 2.5e10，两档与 demand 编译器共用）
#   CCDG_DEMAND_OPTS  ccdg_demand.py 附加参数（如 "--phase" 启用相位格点档）
#   CCDG_INJECT_QDEPTH / CCDG_COMPUTE_RATE 透传 run_ccdg_mesh.sh
#   WSE_PLAN_CAPTURE  LAMMPS 源码级 CommBrick plan 导出（默认 1；0=禁用）
#
# 流程（真实 LAMMPS 只跑一次）:
#   ① 生成 in.lammps（无 minimize，仅 run 段）
#   ② mpirun + DUMPI(LD_PRELOAD) 捕获
#   ③ dumpi2ccdg ×3: trace(原始) / compact(裁剪+折叠) / trimonly(裁剪不折叠)
#   ④ 质量闸门（替代 SimGrid）: BARRIER 锚定全覆盖 + comm_bytes 守恒 +
#      每.rank 首节点 == BARRIER + Loop time 存在
#   ⑤ 按模式注入 BookSim 方形 2D mesh（run_ccdg_mesh.sh 零改动复用）
#   ⑥ evaluation: 解析 stats/log → evaluation.txt + evaluation.json
#      （账本: ranks×makespan = compute + blocked + congestion + sched_wait + idle）
#
# 退出码: 0=全部档 PASS   1=有档 unresolved≠0 或 sent≠recv   2=流水线错误
# 产物: runs/pipeline/<mode>_<体系>_<原子数>a_<rank数>r_<时间戳>/
#   ├── in.lammps / lammps.log / dumpi-*.bin|meta
#   ├── trace_<R>ranks_global.ccdg / compact_<R>ranks_global.ccdg
#   ├── trimonly_<R>ranks_global.ccdg            （两档共同载体，同源对比的前提）
#   ├── quality_gate.txt                          （闸门证据，BARRIER 锚定 span 等）
#   ├── booksim_free/  + booksim_free_result.txt  [free|both]
#   ├── booksim_demand/ + demand_plan.log + booksim_demand_result.txt [demand|both]
#   └── evaluation.txt / evaluation.json
# ============================================================================
set -uo pipefail

ROOT="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
INSTALL=$ROOT/install
LMP=$INSTALL/bin/lmp
LIBDUMPI=$INSTALL/lib/libdumpi.so
DUMPI2CCDG=$ROOT/dumpi2ccdg/dumpi2ccdg
DEMAND_CC=$ROOT/booksim2/ccdg_demand.py
BS_RUNNER=$ROOT/booksim2/run_ccdg_mesh.sh
RESULTS_DIR=$ROOT/booksim2/results
WSE_PLAN_MERGER=$ROOT/merge_wse_plan.py
WSE_PLAN_VALIDATOR=$ROOT/validate_wse_plan.py

MODE=${1:-}
RANKS=${2:-}
SYSTEM=${3:-}
NATOMS=${4:-}
SIMMODE=${5:-both}
CAPTURE_STEPS=${CAPTURE_STEPS:-1}
BS_TIMEOUT=${BOOKSIM_TIMEOUT:-5400}
COMPUTE_CAP=${CCDG_COMPUTE_CAP:-2.5e10}
DEMAND_OPTS=${CCDG_DEMAND_OPTS:-}
WSE_PLAN_CAPTURE=${WSE_PLAN_CAPTURE:-1}

usage() { sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }
log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# ── 参数与组件校验 ────────────────────────────────────────────────────────
case "$MODE" in short|long) ;; *) usage ;; esac
[[ "$RANKS" =~ ^[0-9]+$ ]] && [ "$RANKS" -ge 1 ] || { echo "ERROR: rank数须为正整数"; exit 2; }
case "$SYSTEM" in cu|h2o|lialocl) ;; *) echo "ERROR: 体系须为 cu|h2o|lialocl"; usage ;; esac
[[ "$NATOMS" =~ ^[0-9]+$ ]] && [ "$NATOMS" -ge 1 ] || { echo "ERROR: 原子数须为正整数"; exit 2; }
case "$SIMMODE" in free|demand|both) ;; *) echo "ERROR: 模式须为 free|demand|both"; usage ;; esac
case "$WSE_PLAN_CAPTURE" in 0|1) ;; *) echo "ERROR: WSE_PLAN_CAPTURE 须为 0 或 1"; exit 2 ;; esac
[[ "$CAPTURE_STEPS" =~ ^[0-9]+$ ]] && [ "$CAPTURE_STEPS" -ge 1 ] || { echo "ERROR: CAPTURE_STEPS 须为正整数"; exit 2; }
K=$(awk -v r="$RANKS" 'BEGIN{k=int(sqrt(r)+0.5); if(k*k==r) print k; else print 0}')
[ "$K" -gt 0 ] || { echo "ERROR: rank数 $RANKS 须为完全平方数（BookSim 方形 mesh k=√N）"; exit 2; }

for f in "$LMP" "$LIBDUMPI" "$DUMPI2CCDG" "$DEMAND_CC" "$BS_RUNNER"; do
  [ -e "$f" ] || { echo "ERROR: 缺少组件 $f"; exit 2; }
done
if [ "$WSE_PLAN_CAPTURE" = 1 ]; then
  for f in "$WSE_PLAN_MERGER" "$WSE_PLAN_VALIDATOR"; do
    [ -f "$f" ] || { echo "ERROR: 缺少 WSE plan 工具 $f"; exit 2; }
  done
fi

PAIR_SUFFIX=cut; [ "$MODE" = long ] && PAIR_SUFFIX=long
RUN_DIR=$ROOT/runs/pipeline/${MODE}_${SYSTEM}_${NATOMS}a_${RANKS}r_$(date +%Y%m%d_%H%M%S)
mkdir -p "$RUN_DIR"
log "mode=$MODE sys=$SYSTEM atoms=$NATOMS ranks=$RANKS mesh=${K}x${K} 捕获步数=$CAPTURE_STEPS 注入模式=$SIMMODE cap=$COMPUTE_CAP"
log "运行目录: $RUN_DIR"

# ── ① 生成 in.lammps（与旧流水线同口径：无 minimize，仅 run 段）──────────
make_input() {
  local steps=$1 outdir=$2
  case "$SYSTEM" in
  cu)
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

# Z 压缩 (Cerebras 范式): 2D 处理器网格, Z 邻居本地化
processors      ${K} ${K} 1

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
    if [ $((actual % 2)) -eq 1 ]; then actual=$((actual-1)); fi
    cut=$(awk -v n="$n" 'BEGIN{c=2*n-0.3; if(c>10)c=10; printf "%.2f",c}')
    cat > "$outdir/in.lammps" <<EOF
# H2O 类 2-型 LJ+Coul 体系: target=${NATOMS}, actual=${actual} (fcc ${n}^3*4), mode=${MODE}
units           metal
boundary        p p p
atom_style      charge

# Z 压缩 (Cerebras 范式): 2D 处理器网格, Z 邻居本地化
processors      ${K} ${K} 1

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

# Z 压缩 (Cerebras 范式): 2D 处理器网格, Z 邻居本地化
processors      ${K} ${K} 1

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

# ── ② DUMPI capture ──────────────────────────────────────────────────────
export LD_PRELOAD=$LIBDUMPI
export LD_LIBRARY_PATH=$INSTALL/lib:${LD_LIBRARY_PATH:-}
export DUMPI_OUTDIR=$RUN_DIR
WSE_MPI_ARGS=()
WSE_PLAN_PREFIX=$RUN_DIR/wse_plan
WSE_PLAN=$RUN_DIR/wse_plan.json
if [ "$WSE_PLAN_CAPTURE" = 1 ]; then
  export LAMMPS_WSE_PLAN=$WSE_PLAN_PREFIX
  WSE_MPI_ARGS=(-x LAMMPS_WSE_PLAN)
else
  unset LAMMPS_WSE_PLAN
fi

log "运行真实 LAMMPS (mpirun -np $RANKS, DUMPI 截获中)..."
( cd "$RUN_DIR" && timeout 1800 mpirun -np "$RANKS" --allow-run-as-root --oversubscribe \
    -x LD_PRELOAD -x LD_LIBRARY_PATH -x DUMPI_OUTDIR "${WSE_MPI_ARGS[@]}" \
    "$LMP" -in in.lammps > lammps.log 2>&1 )
RC=$?
if [ $RC -ne 0 ]; then
  echo "ERROR: LAMMPS 退出码 $RC，日志尾部:" >&2
  tail -20 "$RUN_DIR/lammps.log" >&2
  exit 2
fi
NMETA=$(ls "$RUN_DIR"/dumpi-*.meta 2>/dev/null | wc -l)
[ "$NMETA" -ge 1 ] || { echo "ERROR: 未捕获到 DUMPI trace"; exit 2; }
LOOP_TIME=$(grep -oP 'Loop time of \K[0-9.]+' "$RUN_DIR/lammps.log" | tail -1)
[ -n "$LOOP_TIME" ] || { echo "ERROR: log 中无 Loop time（裁剪窗口无法锚定）"; exit 2; }
log "DUMPI meta=$NMETA  Loop time=${LOOP_TIME}s (run 段，纯迭代口径)"
if [ "$WSE_PLAN_CAPTURE" = 1 ]; then
  python3 "$WSE_PLAN_MERGER" "$WSE_PLAN_PREFIX" \
    --expected-ranks "$RANKS" -o "$WSE_PLAN" \
    > "$RUN_DIR/wse_plan_merge.log" 2>&1 || {
      echo "ERROR: WSE plan 合并失败" >&2
      cat "$RUN_DIR/wse_plan_merge.log" >&2
      exit 2
    }
  log "$(cat "$RUN_DIR/wse_plan_merge.log")"
fi

# ── ③ dumpi2ccdg ×3 ─────────────────────────────────────────────────────
CCDG=$RUN_DIR/trace_${RANKS}ranks_global.ccdg
CCDG_COMPACT=$RUN_DIR/compact_${RANKS}ranks_global.ccdg
CCDG_TRIMONLY=$RUN_DIR/trimonly_${RANKS}ranks_global.ccdg

log "生成原始 CCDG ..."
LD_LIBRARY_PATH=$INSTALL/lib:${LD_LIBRARY_PATH:-} \
  "$DUMPI2CCDG" "$RUN_DIR" > "$CCDG" 2> "$RUN_DIR/ccdg_gen.log"
grep -q '"num_ranks"' "$CCDG" || { echo "ERROR: CCDG 生成失败"; tail -5 "$RUN_DIR/ccdg_gen.log"; exit 2; }

log "生成 compact CCDG (CCDG_TRIM_SETUP + CCDG_COMPACT) ..."
CCDG_TRIM_SETUP=1 CCDG_COMPACT=1 LD_LIBRARY_PATH=$INSTALL/lib:${LD_LIBRARY_PATH:-} \
  "$DUMPI2CCDG" "$RUN_DIR" > "$CCDG_COMPACT" 2> "$RUN_DIR/ccdg_compact_gen.log"
grep -q '"num_ranks"' "$CCDG_COMPACT" || { echo "ERROR: compact CCDG 生成失败"; tail -5 "$RUN_DIR/ccdg_compact_gen.log"; exit 2; }

log "生成 trimonly CCDG (CCDG_TRIM_SETUP, 无 COMPACT) ..."
env -u CCDG_COMPACT CCDG_TRIM_SETUP=1 LD_LIBRARY_PATH=$INSTALL/lib:${LD_LIBRARY_PATH:-} \
  "$DUMPI2CCDG" "$RUN_DIR" > "$CCDG_TRIMONLY" 2> "$RUN_DIR/ccdg_trimonly_gen.log"
grep -q '"num_ranks"' "$CCDG_TRIMONLY" || { echo "ERROR: trimonly CCDG 生成失败"; tail -5 "$RUN_DIR/ccdg_trimonly_gen.log"; exit 2; }
log "CCDG: raw=$(wc -c < "$CCDG")B compact=$(wc -c < "$CCDG_COMPACT")B trimonly=$(wc -c < "$CCDG_TRIMONLY")B"

# ── ④ 质量闸门（替代 SimGrid）───────────────────────────────────────────
# 证据 1: 每 rank 的迭代窗口都锚定在 run 首尾 BARRIER 对上（v2 口径）
# 证据 2: compact 折叠 comm_bytes 精确守恒
# 证据 3: 裁剪后每 rank 首节点为 BARRIER（循环起始），无 setup BCAST 残留
# 证据 4: 锚定 span 与 log Loop time 偏差（rank 间 barrier 入口 skew，信息项不判 FAIL）
{
  echo "# 质量闸门 $(date '+%F %T')  载体: trimonly_${RANKS}ranks_global.ccdg"
  echo "loop_time_sec = $LOOP_TIME"
} > "$RUN_DIR/quality_gate.txt"

ANCHORED=$(grep -c "step anchored on BARRIER pair" "$RUN_DIR/ccdg_trimonly_gen.log" || true)
SPAN_DEV=$(grep -oP 'span \K[0-9.]+' "$RUN_DIR/ccdg_trimonly_gen.log" \
           | awk -v lt="$LOOP_TIME" '{d=($1*1e-6-lt); if(d<0)d=-d; s+=d} END{if(NR)printf "%.1f",s/NR*100/lt}')
CONSERVED=$(grep -c "comm bytes conserved" "$RUN_DIR/ccdg_compact_gen.log" || true)
GATE_JSON=$(python3 - "$CCDG_TRIMONLY" <<'PYEOF'
import json, sys, collections
g = json.load(open(sys.argv[1]))
by = collections.defaultdict(list)
for nd in g["nodes"]:
    by[nd["rank"]].append(nd)
first = collections.Counter(v[0]["type"] for v in by.values())
types = collections.Counter(nd["type"] for nd in g["nodes"])
sends = [sum(1 for nd in v if nd["type"] in ("SEND", "ISEND")) for v in by.values()]
import json as j
print(j.dumps({
    "num_ranks": g["num_ranks"],
    "nodes": len(g["nodes"]),
    "cross_edges": len(g.get("cross_rank_edges", [])),
    "first_node_types": dict(first),
    "all_first_barrier": set(first) == {"BARRIER"},
    "setup_bcast_left": types.get("BCAST", 0),
    "send_per_rank_min": min(sends), "send_per_rank_max": max(sends),
    "type_counts": dict(types),
}))
PYEOF
)
ALL_ANCHORED=$([ "$ANCHORED" -eq "$RANKS" ] && echo true || echo false)
BYTES_OK=$([ "$CONSERVED" -ge 1 ] && echo true || echo false)
echo "barrier_anchored_ranks = $ANCHORED / $RANKS"       >> "$RUN_DIR/quality_gate.txt"
echo "anchor_span_vs_loop_avg_pct = ${SPAN_DEV:-NA}%"    >> "$RUN_DIR/quality_gate.txt"
echo "comm_bytes_conserved = $BYTES_OK"                  >> "$RUN_DIR/quality_gate.txt"
echo "$GATE_JSON"                                        >> "$RUN_DIR/quality_gate.txt"

if [ "$WSE_PLAN_CAPTURE" = 1 ]; then
  if [ "$MODE" = long ]; then
    python3 "$WSE_PLAN_VALIDATOR" "$WSE_PLAN" "$CCDG_TRIMONLY" \
      --threshold 0.02 > "$RUN_DIR/wse_plan_validation.json" || true
    echo "wse_plan_vs_ccdg = PARTIAL (CommBrick only; Kspace/FFT not exported)" \
      >> "$RUN_DIR/quality_gate.txt"
    log "WSE plan 对拍: PARTIAL（long 的 Kspace/FFT 尚未纳入 Phase 0 导出）"
  elif python3 "$WSE_PLAN_VALIDATOR" "$WSE_PLAN" "$CCDG_TRIMONLY" \
        --threshold 0.02 > "$RUN_DIR/wse_plan_validation.json"; then
    echo "wse_plan_vs_ccdg = PASS (per-direction bytes <= 2%)" \
      >> "$RUN_DIR/quality_gate.txt"
    log "WSE plan 对拍: PASS (逐方向字节偏差 <= 2%)"
  else
    echo "wse_plan_vs_ccdg = FAIL" >> "$RUN_DIR/quality_gate.txt"
    echo "ERROR: WSE plan 与 trimonly CCDG 对拍失败: $RUN_DIR/wse_plan_validation.json" >&2
    exit 2
  fi
fi

log "闸门: BARRIER锚定 $ANCHORED/$RANKS  span偏差均值 ${SPAN_DEV:-NA}%  bytes守恒=$BYTES_OK"
grep -q '"all_first_barrier": true' <<<"$GATE_JSON" || log "警告: 存在首节点非 BARRIER 的 rank（窗口锚定可疑），详见 quality_gate.txt"
if [ "$ALL_ANCHORED" != true ] || [ "$BYTES_OK" != true ]; then
  grep -E "no boundary|degenerate|rolled back" "$RUN_DIR/ccdg_trimonly_gen.log" | head -5 | sed 's/^/  /' >&2
  log "警告: 闸门证据不完整，结果可用性自担（建议检查 trace 捕获环境）"
fi

# Z 压缩闸门: 处理器网格必须 = K×K×1；short 模式方向集必须 b≤2（z 投影长程类不得残留）
GRID_LINE=$(grep -E "^[[:space:]]*[0-9]+ by [0-9]+ by [0-9]+ MPI processor grid" "$RUN_DIR/log.lammps" | tail -1 | sed 's/^[[:space:]]*//')
echo "proc_grid = ${GRID_LINE:-NONE}" >> "$RUN_DIR/quality_gate.txt"
if [ "$GRID_LINE" != "$K by $K by 1 MPI processor grid" ]; then
  echo "ERROR: 处理器网格非 ${K}x${K}x1（Z 压缩口径破坏）: '${GRID_LINE:-NONE}'" >&2
  exit 2
fi
if [ "$MODE" = short ]; then
  DIRCHK=$(python3 - "$CCDG_TRIMONLY" <<'PYEOF'
import json, sys, collections
g = json.load(open(sys.argv[1]))
k = int(round(g["num_ranks"] ** 0.5))
allowed = {(1, 0), (-1, 0), (0, 1), (0, -1),
           (k - 1, 0), (-(k - 1), 0), (0, k - 1), (0, -(k - 1))}
dirs = collections.Counter()
for nd in g["nodes"]:
    if nd["type"] in ("SEND", "ISEND"):
        d = nd.get("comm_dst", nd["rank"])
        dirs[(d % k - nd["rank"] % k, d // k - nd["rank"] // k)] += 1
bad = {d: c for d, c in dirs.items() if d not in allowed}
print(f"{len(dirs)} {len(bad)} " + repr(dict(sorted(bad.items()))))
PYEOF
)
  NDIR=$(echo "$DIRCHK" | awk '{print $1}')
  NBAD=$(echo "$DIRCHK" | awk '{print $2}')
  echo "dir_classes = $NDIR  dir_off_lattice = $NBAD  (须为面邻居±1或PBC接缝±(k-1))" >> "$RUN_DIR/quality_gate.txt"
  [ "$NBAD" -eq 0 ] || { echo "ERROR: 方向集含 2D 面邻居/接缝之外的类（z 投影残留），Z 压缩失败: $DIRCHK" >&2; exit 2; }
  log "Z 压缩闸门: grid=${K}x${K}x1  dir_classes=$NDIR (全部面邻居/接缝)  PASS"
fi

# ── ⑤ BookSim 注入（run_ccdg_mesh.sh 零改动复用；TAG=父目录名，两档隔离）──
declare -a RAN_TAGS=()
inject_free() {
  local d=$RUN_DIR/booksim_free
  mkdir -p "$d"
  cp "$CCDG_TRIMONLY" "$d/trimonly_${RANKS}ranks_global.ccdg"
  log "BookSim free 档 (trimonly 同源载体, 1 步注入) ..."
  CCDG_COMPUTE_CAP=$COMPUTE_CAP "$BS_RUNNER" "$d/trimonly_${RANKS}ranks_global.ccdg" \
      "$BS_TIMEOUT" 1 2>&1 | tee "$RUN_DIR/booksim_free_result.txt"
  RAN_TAGS+=("free:$RUN_DIR/booksim_free_result.txt")
}

inject_demand() {
  local d=$RUN_DIR/booksim_demand
  mkdir -p "$d"
  log "需求提炼 + 时槽排布 + 无同步 CCDG (ccdg_demand.py ${DEMAND_OPTS:-}) ..."
  python3 "$DEMAND_CC" "$CCDG_TRIMONLY" \
    -o "$d/demand_${RANKS}ranks" \
    --cap "$COMPUTE_CAP" $DEMAND_OPTS 2>&1 | tee "$RUN_DIR/demand_plan.log" || exit 2
  DEMAND_CCDG=$d/demand_${RANKS}ranks_demand.ccdg
  DEMAND_EST=$d/demand_${RANKS}ranks_demand.est
  grep -q "plan report" "$RUN_DIR/demand_plan.log" && [ -f "$DEMAND_CCDG" ] && [ -f "$DEMAND_EST" ] || {
    echo "ERROR: demand CCDG 产物缺失"; exit 2; }
  log "BookSim demand 档 (无同步 CCDG + est 门控, 1 步注入) ..."
  CCDG_COMPUTE_CAP=$COMPUTE_CAP CCDG_SCHED_FILE="$DEMAND_EST" \
    "$BS_RUNNER" "$DEMAND_CCDG" "$BS_TIMEOUT" 1 \
    2>&1 | tee "$RUN_DIR/booksim_demand_result.txt"
  RAN_TAGS+=("demand:$RUN_DIR/booksim_demand_result.txt")
}

declare -a RAN_TAGS=()
inject_free() {
  local d=$RUN_DIR/booksim_free
  mkdir -p "$d"
  cp "$CCDG_TRIMONLY" "$d/trimonly_${RANKS}ranks_global.ccdg"
  log "BookSim free 档 (logical 载体 free 注入, 1 步) ..."
  CCDG_COMPUTE_CAP=$COMPUTE_CAP "$BS_RUNNER" "$d/trimonly_${RANKS}ranks_global.ccdg" \
      "$BS_TIMEOUT" 1 2>&1 | tee "$RUN_DIR/booksim_free_result.txt"
  RAN_TAGS+=("free:$RUN_DIR/booksim_free_result.txt")
}

inject_free_fold() {
  local d=$RUN_DIR/booksim_free_fold
  mkdir -p "$d"
  log "free_fold: PBC 交叉排布重标号 (paper III-E fold) ..."
  python3 "$DEMAND_CC" --fold-only "$CCDG_TRIMONLY" \
    "$d/trimonly_${RANKS}ranks_global.ccdg" || exit 2
  log "BookSim free_fold 档 (fold 载体 free 注入, 1 步) ..."
  CCDG_COMPUTE_CAP=$COMPUTE_CAP "$BS_RUNNER" \
    "$d/trimonly_${RANKS}ranks_global.ccdg" "$BS_TIMEOUT" 1 \
    2>&1 | tee "$RUN_DIR/booksim_free_fold_result.txt"
  RAN_TAGS+=("free_fold:$RUN_DIR/booksim_free_fold_result.txt")
}

# compile_and_inject <tag> <ccdg_demand.py 附加参数...>
compile_and_inject() {
  local tag=$1; shift
  local d=$RUN_DIR/booksim_$tag
  mkdir -p "$d"
  log "demand[$tag] 编译 (ccdg_demand.py $*) ..."
  python3 "$DEMAND_CC" "$CCDG_TRIMONLY" \
    -o "$d/demand_${RANKS}ranks" \
    --cap "$COMPUTE_CAP" "$@" 2>&1 | tee "$RUN_DIR/demand_${tag}_plan.log" || exit 2
  local ccdg=$d/demand_${RANKS}ranks_demand.ccdg
  local est=$d/demand_${RANKS}ranks_demand.est
  grep -q "plan report" "$RUN_DIR/demand_${tag}_plan.log" \
    && [ -f "$ccdg" ] && [ -f "$est" ] || {
    echo "ERROR: demand[$tag] 产物缺失"; exit 2; }
  log "BookSim demand[$tag] 档 (est 门控, 1 步注入) ..."
  CCDG_COMPUTE_CAP=$COMPUTE_CAP CCDG_SCHED_FILE="$est" \
    "$BS_RUNNER" "$ccdg" "$BS_TIMEOUT" 1 \
    2>&1 | tee "$RUN_DIR/booksim_${tag}_result.txt"
  RAN_TAGS+=("$tag:$RUN_DIR/booksim_${tag}_result.txt")
}

case "$SIMMODE" in
  free)   inject_free ;;
  demand) compile_and_inject hb $DEMAND_OPTS --fold-pbc
          compile_and_inject cerebras $DEMAND_OPTS --fold-pbc --cerebras
          compile_and_inject ilv $DEMAND_OPTS --fold-pbc --cerebras --wavelet-mode interleave ;;
  both)   inject_free
          inject_free_fold
          compile_and_inject hb $DEMAND_OPTS --fold-pbc
          compile_and_inject cerebras $DEMAND_OPTS --fold-pbc --cerebras
          compile_and_inject ilv $DEMAND_OPTS --fold-pbc --cerebras --wavelet-mode interleave ;;
esac

# ── ⑥ evaluation：stats/log → evaluation.txt + evaluation.json ───────────
python3 - "$RUN_DIR" "$RANKS" "$CAPTURE_STEPS" "$LOOP_TIME" "$SIMMODE" "$COMPUTE_CAP" \
         "${RAN_TAGS[@]:-}" <<'PYEOF'
import json, os, re, sys

run_dir, ranks, cap_steps, loop_s, simmode, cap = sys.argv[1:7]
tags = sys.argv[7:]
ranks = int(ranks); cap_steps = int(cap_steps); loop_s = float(loop_s)
NOC_NS = 0.5  # 2.0 GHz

def parse_result(txt):
    """从 run_ccdg_mesh.sh 的 teed 输出定位 cfg/stats/log 并提取评估字段"""
    if not txt or not os.path.exists(txt):
        return None
    body = open(txt, errors="replace").read()
    m = re.search(r"cfg=(\S+)", body)
    if not m:
        return None
    cfg = m.group(1)
    stem = cfg[:-len(".cfg")] if cfg.endswith(".cfg") else cfg
    stats_p, log_p = stem + "_stats.txt", stem + ".log"
    d = {"cfg": os.path.basename(cfg)}
    for f in (stats_p, log_p):
        if os.path.exists(f):
            d[os.path.basename(f)] = open(f, errors="replace").read()
    stats = d.get(os.path.basename(stats_p), "")
    log_b = d.get(os.path.basename(log_p), "")
    def g(pat, src, cast=float):
        mm = re.findall(pat, src)
        return cast(mm[-1]) if mm else None
    out = {
        "tag": os.path.basename(stem),
        "total_cycles": g(r"total_sim_cycles = (\d+)", stats, int)
                        or g(r"completed in (\d+)", log_b, int),
        "compute_cycles": g(r"compute_cycles = (\d+)", stats, int) or 0,
        "blocked_cycles": g(r"blocked_cycles = (\d+)", stats, int) or 0,
        "congestion_cycles": g(r"congestion_cycles = (\d+)", stats, int) or 0,
        "sched_wait_cycles": g(r"sched_wait_cycles = (\d+)", stats, int) or 0,
        "packets_sent": g(r"Packets sent: (\d+)", log_b, int),
        "packets_recv": g(r"received: (\d+)", log_b, int),
        "unresolved": g(r"WARNING: (\d+) cross-rank edges", log_b, int) or 0,
    }
    if out["total_cycles"] is None:
        return None
    return out

sims = {}
missing = []
for t in tags:
    name, path = t.split(":", 1)
    r = parse_result(path)
    if r:
        sims[name] = r
    else:
        missing.append(name)
        sims[name] = {"tag": name, "total_cycles": 0, "compute_cycles": 0,
                      "blocked_cycles": 0, "congestion_cycles": 0, "sched_wait_cycles": 0,
                      "packets_sent": None, "packets_recv": None, "unresolved": -1,
                      "wall_sec_at_2ghz": 0.0, "cycles_per_iter": 0, "timesteps_per_sec": 0,
                      "ledger_per_rank_avg": {k: 0 for k in ("compute", "blocked", "congestion", "sched_wait")},
                      "ledger_ratio": {k: 0.0 for k in ("compute", "blocked", "congestion", "sched_wait")},
                      "pe_done_idle_ratio": 0.0,
                      "verdict": "FAIL", "error": f"未能从 {path} 解析出 BookSim 结果"}
        continue

for name, s in sims.items():
    c = s["total_cycles"]
    if not c:
        continue
    s["wall_sec_at_2ghz"] = c * NOC_NS / 1e9
    s["cycles_per_iter"] = round(c / cap_steps)
    s["timesteps_per_sec"] = round(1e9 / (s["cycles_per_iter"] * NOC_NS))
    mk = ranks * c  # 账本: ranks×makespan（PE-cycle 总量）
    led = {k: s[f"{k}_cycles"] for k in ("compute", "blocked", "congestion", "sched_wait")}
    s["ledger_per_rank_avg"] = {k: round(v / ranks) for k, v in led.items()}
    s["ledger_ratio"] = {k: round(v / mk, 4) for k, v in led.items()}
    s["pe_done_idle_ratio"] = round(1.0 - sum(led.values()) / mk, 4)
    s["verdict"] = ("PASS" if s["unresolved"] == 0
                    and s["packets_sent"] in (None, s["packets_recv"])
                    else "FAIL")

comparisons = []
PAIR_DEFS = [("free_fold", "hb", "需求驱动净收益(正确性口径, fold 载体)"),
             ("hb", "cerebras", "XY stage 串行化+barrier 代价(PBC fold 后)"),
             ("cerebras", "ilv", "相位串行 vs 交错代价"),
             ("free", "free_fold", "PBC 交叉排布代价(反应式档)")]
for a, b, note in PAIR_DEFS:
    if (a in sims and b in sims and sims[a].get("total_cycles")
            and sims[b].get("total_cycles")):
        d = (sims[a]["total_cycles"] - sims[b]["total_cycles"]) \
            / sims[a]["total_cycles"] * 100
        comparisons.append({"from": a, "to": b,
                            "a_cycles": sims[a]["total_cycles"],
                            "b_cycles": sims[b]["total_cycles"],
                            "delta_pct": round(d, 2), "note": note})

result = {
    "run_dir": run_dir, "ranks": ranks, "capture_steps": cap_steps,
    "noc_ghz": 2.0, "compute_cap_ops": cap, "real_loop_time_sec": loop_s,
    "carrier": "trimonly", "simmode": simmode, "sims": sims,
    "comparisons": comparisons,
}
json.dump(result, open(os.path.join(run_dir, "evaluation.json"), "w"),
          indent=2, ensure_ascii=False)

# ── 人类可读报告 ──
L = []
L.append(f"==== NoC pipeline evaluation ====")
L.append(f"run_dir={run_dir}")
L.append(f"ranks={ranks} capture_steps={cap_steps} noc=2.0GHz cap={cap} ops/s")
L.append(f"real_loop_time={loop_s*1e6:.1f} us (LAMMPS run 段, 真实集群口径, 仅作参照)")
hdr = f"{'档':<8}{'cycles':>12}{'cyc/iter':>10}{'steps/s':>9}{'blocked%':>9}{'cong%':>7}{'sched%':>7}{'idle%':>7}{'unres':>6}{'verdict':>8}"
for name, s in sims.items():
    lr = s["ledger_ratio"]
    L.append(hdr)
    L.append(f"{name:<8}{s['total_cycles']:>12,}{s['cycles_per_iter']:>10,}"
             f"{s['timesteps_per_sec']:>9,}{lr['blocked']*100:>8.2f}%{lr['congestion']*100:>6.2f}%"
             f"{lr['sched_wait']*100:>6.2f}%{s['pe_done_idle_ratio']*100:>6.2f}%"
             f"{s['unresolved']:>6}{s['verdict']:>8}")
    L.append(f"        账本(rank均值): compute={s['ledger_per_rank_avg']['compute']:,} "
             f"blocked={s['ledger_per_rank_avg']['blocked']:,} "
             f"congestion={s['ledger_per_rank_avg']['congestion']:,} "
             f"sched_wait={s['ledger_per_rank_avg']['sched_wait']:,} "
             f"sent={s['packets_sent']} recv={s['packets_recv']}")
    L.append(f"        {s['tag']}")
if comparisons:
    L.append(f"---- 同源档位对比 ----")
    for cp in comparisons:
        L.append(f"{cp['from']} -> {cp['to']}: "
                 f"{cp['a_cycles']:,} -> {cp['b_cycles']:,} "
                 f"({cp['delta_pct']:+.2f}%)  [{cp['note']}]")
overall = all(s["verdict"] == "PASS" for s in sims.values()) and sims
L.append(f"OVERALL: {'PASS' if overall else 'FAIL'}")
open(os.path.join(run_dir, "evaluation.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
sys.exit(0 if overall else 1)
PYEOF
EVAL_RC=$?

log "完成。评估: $RUN_DIR/evaluation.txt (+.json)"
exit $EVAL_RC
