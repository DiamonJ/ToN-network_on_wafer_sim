#!/bin/bash
# WSE wavelet 主实验：3 CCDG × {free, phase_width 1/16/64/256}（strip_width=2, num_vcs=24）
set -uo pipefail
BS=/work1/jiangtao/lammps_trace/booksim2
R16_SHORT=/work1/jiangtao/lammps_trace/runs/pipeline/short_lialocl_2688a_16r_20260821_163455/compact_16ranks_global.ccdg
R64_SHORT=/work1/jiangtao/lammps_trace/runs/pipeline/short_lialocl_2688a_64r_20260821_171510/compact_64ranks_global.ccdg
R16_LONG=/work1/jiangtao/lammps_trace/runs/pipeline/long_lialocl_2688a_16r_20260821_163513/compact_16ranks_global.ccdg

run() {
    local ccdg=$1 gating=$2 pw=$3
    local name; name="$(basename "$(dirname "$ccdg")")"
    echo "=== $(date +%H:%M:%S) $name gating=$gating phase_width=$pw ==="
    bash "$BS/run_ccdg_mesh.sh" "$ccdg" 1200 1 "$gating" "$pw" 2 || echo "FAIL rc=$?"
}

# free 基线已在 run_ccdg_mesh_results.csv 中验证（本次会话回归确认一致），
# 这里只重跑 wse 模式（修复死锁后的有效结果）
for pw in 1 16 64 256; do
    run "$R16_SHORT" 1 "$pw"
done

for pw in 1 16 64 256; do
    run "$R64_SHORT" 1 "$pw"
done

for pw in 1 16 64 256; do
    run "$R16_LONG" 1 "$pw"
done

echo "=== ALL DONE $(date +%H:%M:%S) ==="
