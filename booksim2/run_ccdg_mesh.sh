#!/bin/bash
# run_ccdg_mesh.sh — 把 CCDG 注入 BookSim 方形 2D mesh 仿真
# 用法: ./run_ccdg_mesh.sh <ccdg文件> [超时秒数,默认3600] [步数,默认1] [wse_gating,默认0] [wse_phase_width,默认64] [wse_strip_width,默认2]
# 要求: CCDG 的 num_ranks 为完全平方数 N，mesh 取 k=sqrt(N), n=2
# 频率口径: noc = 2.0 GHz，折算墙钟秒 = cycles × 0.5ns；cpu_frequency_ghz 保留模板 trace 机频率 2.5（仅 legacy 回退与信息展示，v3 计算由 capability 驱动）
# v3 口径: 计算能力用环境变量 CCDG_COMPUTE_CAP 覆盖（默认 2.5e9 ops/s = 2.5GHz×1op/cycle，
#          等价于旧 ccdg_compute_rate=1.25 ops/cycle @ noc 2.0GHz）；算力扫描传 1.25e9(=0.5×)/2.5e9(=1×)/5e9(=2×) 等；
#          也可用 CCDG_COMPUTE_RATE（ops/cycle）直接指定（capability 置 0）
# CSV: ccdg_dir,...,unresolved,steps,cycles_per_iter,timesteps_per_sec,gating,inject_blocked_cycles,blocked_cycles,blocked_ratio
#      (cycles_per_iter=total_cycles/步数; timesteps_per_sec=1e9/(cycles_per_iter×NOC_PERIOD_NS))
set -euo pipefail

BS_DIR="$(cd "$(dirname "$0")" && pwd)"
BOOKSIM="$BS_DIR/booksim"
TEMPLATE="$BS_DIR/ccdg_lammps_4x4.cfg"
RESULTS_DIR="$BS_DIR/results"
CSV="$RESULTS_DIR/ccdg_mesh_results.csv"
NOC_PERIOD_NS=0.5   # 2.0 GHz

CCDG="${1:?用法: $0 <ccdg文件> [超时秒数] [步数] [wse_gating] [wse_phase_width] [wse_strip_width]}"
TIMEOUT="${2:-3600}"
STEPS="${3:-1}"
WSE_GATING="${4:-0}"
WSE_PHASE_WIDTH="${5:-64}"
WSE_STRIP_WIDTH="${6:-2}"
CCDG="$(cd "$(dirname "$CCDG")" && pwd)/$(basename "$CCDG")"
[ -f "$CCDG" ] || { echo "错误: CCDG 文件不存在: $CCDG" >&2; exit 2; }
[ -x "$BOOKSIM" ] || { echo "错误: booksim 未编译: $BOOKSIM" >&2; exit 2; }

mkdir -p "$RESULTS_DIR"

# 1) 读 num_ranks，校验完全平方
PARSE_OUT="$(python3 - "$CCDG" <<'EOF'
import json, math, sys
d = json.load(open(sys.argv[1]))
n = d["num_ranks"]
k = int(round(math.sqrt(n)))
if k * k != n:
    print(f"错误: num_ranks={n} 不是完全平方数，BookSim mesh(k^n) 无法构造方形 2D", file=sys.stderr)
    sys.exit(1)
print(n, k)
EOF
)" || exit 2
read -r RANKS K <<< "$PARSE_OUT"

# 2) 从运行目录名推断 mode（short/long）
NAME="$(basename "$CCDG" .ccdg)"            # trace_<r>ranks_global
DIRTAG="$(basename "$(dirname "$CCDG")")"   # {mode}_{sys}_{a}a_{r}r_{ts}
MODE="${DIRTAG%%_*}"
MESH="${K}x${K}"
[ "$WSE_GATING" = "0" ] && GATING_TAG="free" || GATING_TAG="wse_p${WSE_PHASE_WIDTH}_w${WSE_STRIP_WIDTH}"
# 计算能力参数（环境变量覆盖，默认 2.5e9 ops/s = 2.5GHz×1op/cycle）
COMPUTE_CAP="${CCDG_COMPUTE_CAP:-2.5e10}"
COMPUTE_RATE="${CCDG_COMPUTE_RATE:-0.0}"  # 注意必须带小数点: 整数 0 会走 Assign(int) 而 ccdg_compute_rate 只在 _float_map
# 注入队列深度（flit 数，PE 反压）；0=auto=最大单包长度；可覆盖为任意正整数
INJ_QDEPTH="${CCDG_INJECT_QDEPTH:-0}"
# 确定性编排（EST 表，可选）：非空时 PE 按编排器给出的 release time 执行
SCHED_FILE="${CCDG_SCHED_FILE:-}"
# 非默认队列深度时在命名中附加标识，避免 auto/unlimited 同名覆盖 cfg/stats/log
QD_TAG=""
if [ "$INJ_QDEPTH" != "0" ]; then
    QD_TAG="_qd${INJ_QDEPTH}"
fi
SCHED_TAG=""
if [ -n "$SCHED_FILE" ]; then
    SCHED_TAG="_sched"
    # 相位化 est 表（ccdg_scheduler.py --pw）在文件名中带 _pw<pw>_w<w>，
    # 提取进 TAG 避免不同相位参数互相覆盖 cfg/log/stats/CSV
    PHASE_TAG="$(basename "$SCHED_FILE" | grep -oP '_pw[0-9.eE+]+_w[0-9]+' | head -1 || true)"
    [ -n "$PHASE_TAG" ] && SCHED_TAG="${SCHED_TAG}${PHASE_TAG}"
fi
# booksim yacc 解析器不支持空字符串值（config.l：只认 "name = 'val';" 语法，
# "//" 注释），无 schedule 时把模板行注释掉，勿留 "= ;" 空值否则 Parse error
if [ -n "$SCHED_FILE" ]; then
    SCHED_SED_LINE="ccdg_schedule_file = $SCHED_FILE;"
else
    SCHED_SED_LINE="//ccdg_schedule_file = ;"
fi
# 非默认算力时在命名/CSV 中附加标识，避免不同算力互相覆盖 cfg/stats/log
CAP_TAG=""
if [ "$COMPUTE_CAP" != "2.5e9" ]; then
    CAP_TAG="_cap${COMPUTE_CAP}"
elif [ "$COMPUTE_RATE" != "0.0" ]; then
    CAP_TAG="_rate${COMPUTE_RATE}"
fi
TAG="${DIRTAG}${CAP_TAG}${QD_TAG}${SCHED_TAG}_${STEPS}s_${GATING_TAG}"
CFG="$RESULTS_DIR/ccdg_mesh_${TAG}.cfg"
STATS="$RESULTS_DIR/ccdg_mesh_${TAG}_stats.txt"
LOG="$RESULTS_DIR/ccdg_mesh_${TAG}.log"

# 3) 由模板生成 cfg：覆盖 ccdg_file / 频率 / flit 大小 / k / stats_out / wse 参数
sed -e "s|^ccdg_file = .*|ccdg_file = $CCDG;|" \
    -e "s|^noc_frequency_ghz = .*|noc_frequency_ghz = 2.0;|" \
    -e "s|^flit_size_bytes = .*|flit_size_bytes = 1;|" \
    -e "s|^ccdg_compute_capability = .*|ccdg_compute_capability = $COMPUTE_CAP;|" \
    -e "s|^ccdg_compute_rate = .*|ccdg_compute_rate = $COMPUTE_RATE;|" \
    -e "s|^ccdg_inject_queue_depth = .*|ccdg_inject_queue_depth = $INJ_QDEPTH;|" \
    -e "s|^ccdg_schedule_file = .*|${SCHED_SED_LINE}|" \
    -e "s|^k = .*|k = $K;|" \
    -e "s|^stats_out = .*|stats_out = $STATS;|" \
    -e "s|^wse_gating = .*|wse_gating = $WSE_GATING;|" \
    -e "s|^wse_phase_width = .*|wse_phase_width = $WSE_PHASE_WIDTH;|" \
    -e "s|^wse_strip_width = .*|wse_strip_width = $WSE_STRIP_WIDTH;|" \
    "$TEMPLATE" > "$CFG"

echo "[$(date +%H:%M:%S)] $DIRTAG: ranks=$RANKS mesh=$MESH mode=$MODE steps=$STEPS gating=$GATING_TAG" >&2
echo "[$(date +%H:%M:%S)] cfg=$CFG" >&2

# 4) 运行 BookSim（超时保护）
# 注意: booksim 成功时返回 255（main.cpp: return result ? -1 : 0），仅 124(超时) 视为失败
RC=0
timeout "$TIMEOUT" "$BOOKSIM" "$CFG" > "$LOG" 2>&1 || RC=$?
if [ $RC -eq 124 ]; then
    echo "错误: 仿真超时 (${TIMEOUT}s): $LOG" >&2
fi
[ $RC -ne 0 ] && [ $RC -ne 255 ] && { echo "错误: booksim 异常退出 rc=$RC: $LOG" >&2; }

# 5) 解析结果
CYCLES="$(grep -oP 'completed in \K[0-9]+' "$LOG" | tail -1 || true)"
SENT="$(grep -oP 'Packets sent: \K[0-9]+' "$LOG" | tail -1 || true)"
RECV="$(grep -oP 'received: \K[0-9]+' "$LOG" | tail -1 || true)"
UNRES="$(grep -oP 'WARNING: \K[0-9]+(?= cross-rank edges)' "$LOG" | tail -1 || true)"
UNRES="${UNRES:-0}"
BLOCKED="$(grep -oP 'wse_inject_blocked_cycles = \K[0-9]+' "$STATS" 2>/dev/null | tail -1 || true)"
BLOCKED="${BLOCKED:-}"
if [ -z "$CYCLES" ]; then
    # 回退到 stats 文件
    CYCLES="$(grep -oP 'total_sim_cycles = \K[0-9]+' "$STATS" 2>/dev/null | tail -1 || true)"
fi
if [ -z "$CYCLES" ]; then
    echo "错误: 未找到完成 cycle 数，仿真可能失败，详见 $LOG" >&2
    exit 1
fi
WALL_S="$(awk -v c="$CYCLES" -v p="$NOC_PERIOD_NS" 'BEGIN{printf "%.6f", c*p/1e9}')"
CPI="$(awk -v c="$CYCLES" -v s="$STEPS" 'BEGIN{printf "%.0f", c/s}')"
TPS="$(awk -v c="$CPI" -v p="$NOC_PERIOD_NS" 'BEGIN{printf "%.0f", 1e9/(c*p)}')"

# 6) 追加 CSV（存量旧格式自动迁移: 先补 steps/cycles_per_iter，再补 timesteps_per_sec，再补 gating/inject_blocked_cycles）
if [ -f "$CSV" ] && ! head -1 "$CSV" | grep -q "steps"; then
    awk -F, 'NR==1{print $0",steps,cycles_per_iter,timesteps_per_sec"; next}
             {cpi=$5/100; printf "%s,100,%.0f,%.0f\n", $0, cpi, 1e9/(cpi*0.5)}' "$CSV" > "$CSV.tmp" && mv "$CSV.tmp" "$CSV"
fi
if [ -f "$CSV" ] && ! head -1 "$CSV" | grep -q "timesteps_per_sec"; then
    awk -F, 'NR==1{print $0",timesteps_per_sec"; next}
             {printf "%s,%.0f\n", $0, 1e9/($11*0.5)}' "$CSV" > "$CSV.tmp" && mv "$CSV.tmp" "$CSV"
fi
if [ -f "$CSV" ] && ! head -1 "$CSV" | grep -q "inject_blocked_cycles"; then
    awk -F, 'NR==1{print $0",gating,inject_blocked_cycles"; next}
             {printf "%s,free,\n", $0}' "$CSV" > "$CSV.tmp" && mv "$CSV.tmp" "$CSV"
fi
if [ -f "$CSV" ] && ! head -1 "$CSV" | grep -q "blocked_cycles"; then
    awk -F, 'NR==1{print $0",blocked_cycles,blocked_ratio"; next}
             {printf "%s,,,\n", $0}' "$CSV" > "$CSV.tmp" && mv "$CSV.tmp" "$CSV"
fi
[ -f "$CSV" ] || echo "ccdg_dir,mode,ranks,mesh,total_cycles,wall_sec,packets_sent,packets_recv,unresolved,steps,cycles_per_iter,timesteps_per_sec,gating,inject_blocked_cycles,blocked_cycles,blocked_ratio" > "$CSV"
# 从 stats 文件解析 v3 统计（blocked_cycles/blocked_ratio 仅 v3 CCDG 有）
COMPUTE_CYCLES="$(grep -oP 'compute_cycles = \K[0-9]+' "$STATS" 2>/dev/null | tail -1 || true)"
BLOCKED_CYCLES="$(grep -oP 'blocked_cycles = \K[0-9]+' "$STATS" 2>/dev/null | tail -1 || true)"
BLOCKED_RATIO="$(grep -oP 'blocked_ratio = \K[0-9.]+' "$STATS" 2>/dev/null | tail -1 || true)"
echo "$DIRTAG$CAP_TAG$SCHED_TAG,$MODE,$RANKS,$MESH,$CYCLES,$WALL_S,${SENT:-},${RECV:-},$UNRES,$STEPS,$CPI,$TPS,$GATING_TAG,${BLOCKED:-},${BLOCKED_CYCLES:-},${BLOCKED_RATIO:-}" >> "$CSV"

echo "[$(date +%H:%M:%S)] 完成: cycles=$CYCLES (≈${WALL_S}s @2GHz) steps=$STEPS cycles/iter=$CPI ($TPS steps/s) gating=$GATING_TAG blocked=${BLOCKED:-?} sent=${SENT:-?} recv=${RECV:-?} unresolved=$UNRES rc=$RC" >&2
if [ "$UNRES" != "0" ]; then
    echo "警告: 存在未解析跨 rank 边: $UNRES" >&2
fi
# booksim 成功返回 255；超时 124 原样传播，其余异常码也传播
[ $RC -eq 255 ] && RC=0
exit $RC
