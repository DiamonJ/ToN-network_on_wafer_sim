#!/bin/bash
# 标准 TrafficPattern 交叉验证：{uniform,bitrev,tornado,transpose} × rate{0.1,0.3} × gating{0,1}
set -uo pipefail
BS=/work1/jiangtao/lammps_trace/booksim2
TEMPLATE="$BS/traffic_wse_4x4.cfg"
OUT="$BS/results/traffic_wse_4x4.csv"

[ -f "$OUT" ] || echo "traffic,injection_rate,wse_gating,avg_latency,avg_throughput,inject_blocked" > "$OUT"

run_one() {
    local traffic=$1 rate=$2 gating=$3
    local tag="${traffic}_r${rate}_g${gating}"
    local cfg="$BS/results/traffic_wse_${tag}.cfg"
    local log="$BS/results/traffic_wse_${tag}.log"
    sed -e "s|^traffic = .*|traffic = $traffic;|" \
        -e "s|^injection_rate = .*|injection_rate = $rate;|" \
        -e "s|^wse_gating = .*|wse_gating = $gating;|" \
        -e "s|^stats_out = .*|stats_out = $BS/results/traffic_wse_${tag}_stats.txt;|" \
        "$TEMPLATE" > "$cfg"
    echo "=== $(date +%H:%M:%S) traffic=$traffic rate=$rate gating=$gating ==="
    timeout 600 "$BS/booksim" "$cfg" > "$log" 2>&1
    local rc=$?
    local cycles lat thr blk
    # Latency-mode output has no total-cycle line; use the packet latency
    # and accepted rate lines (sample-period stats, not absolute counters)
    lat="$(grep -oP 'Packet latency average = \K[0-9.]+' "$log" | tail -1 || true)"
    thr="$(grep -oP 'Accepted packet rate average = \K[0-9.]+' "$log" | tail -1 || true)"
    blk="$(grep -oP 'wse_inject_blocked_cycles = \K[0-9]+' "$BS/results/traffic_wse_${tag}_stats.txt" 2>/dev/null | tail -1 || true)"
    echo "$traffic,$rate,$gating,${lat:-},${thr:-},${blk:-}" >> "$OUT"
    echo "lat=${lat:-?} thr=${thr:-?} blocked=${blk:-0} rc=$rc"
}

for traffic in uniform bitrev tornado transpose; do
    for rate in 0.1 0.3; do
        run_one "$traffic" "$rate" 0
        run_one "$traffic" "$rate" 1
    done
done
echo "=== TRAFFIC DONE $(date +%H:%M:%S) ==="
