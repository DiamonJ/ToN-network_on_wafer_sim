#!/bin/bash
# cap_scan_free.sh — wse_gating=0 下扫描算力（ccdg_compute_capability），
# 从 1× (2.5e9 ops/s) 倍增，直到计算耗时占比 < 30% (blocked_ratio > 0.70) 即停。
# 用法: bash cap_scan_free.sh [ccdg文件...] （缺省扫 4 个 v3 CCDG）
set -uo pipefail

BS_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$BS_DIR"

if [ $# -gt 0 ]; then
  CCDGS=("$@")
else
  P=/work1/jiangtao/lammps_trace/runs/pipeline
  CCDGS=(
    "$P/short_lialocl_2688a_4r_20260821_144747/compact_v3.ccdg"
    "$P/short_lialocl_2688a_16r_20260821_163455/compact_v3.ccdg"
    "$P/short_lialocl_2688a_64r_20260821_171510/compact_v3.ccdg"
    "$P/long_lialocl_2688a_16r_20260821_163513/compact_v3.ccdg"
  )
fi

CAPS=(2.5e9 5e9 1e10 2e10 4e10 8e10)   # 1× 2× 4× 8× 16× 32×
THRESH=0.70                             # blocked_ratio 阈值（计算占比 < 30%）

echo "===== 算力扫描 (free gating, 目标 blocked_ratio > $THRESH 即计算占比 < 30%) ====="
for ccdg in "${CCDGS[@]}"; do
  DIRTAG="$(basename "$(dirname "$ccdg")")"
  echo "===== $DIRTAG ====="
  for cap in "${CAPS[@]}"; do
    echo "--- [$(date +%H:%M:%S)] cap=$cap ($(awk -v c=$cap -v b=2.5e9 'BEGIN{printf "%.1fx", c/b}')) ---"
    CCDG_COMPUTE_CAP="$cap" ./run_ccdg_mesh.sh "$ccdg" 600 1 > /dev/null 2>&1
    # 解析该 run 的 stats（cap=2.5e9 时无 _cap 后缀）
    if [ "$cap" = "2.5e9" ]; then
      STATS="$BS_DIR/results/ccdg_mesh_${DIRTAG}_1s_free_stats.txt"
    else
      STATS="$BS_DIR/results/ccdg_mesh_${DIRTAG}_cap${cap}_1s_free_stats.txt"
    fi
    if [ -z "$STATS" ]; then
      echo "!! 未找到 stats，跳过"
      continue
    fi
    BR="$(grep -oP 'blocked_ratio = \K[0-9.]+' "$STATS" | tail -1)"
    TOT="$(grep -oP 'total_sim_cycles = \K[0-9]+' "$STATS" | tail -1)"
    CMP="$(grep -oP 'compute_cycles = \K[0-9]+' "$STATS" | tail -1)"
    BLK="$(grep -oP 'blocked_cycles = \K[0-9]+' "$STATS" | tail -1)"
    echo "  cap=$cap total=$TOT compute=$CMP blocked=$BLK blocked_ratio=$BR 计算占比=$(awk -v r="$BR" 'BEGIN{printf "%.1f%%", (1-r)*100}')"
    OK="$(awk -v r="$BR" -v t="$THRESH" 'BEGIN{print (r>t)?1:0}')"
    if [ "$OK" = "1" ]; then
      echo "  >>> 达标 (计算占比 < 30%)，停止本 CCDG 扫描"
      break
    fi
  done
done
echo "===== 扫描完成 ====="
