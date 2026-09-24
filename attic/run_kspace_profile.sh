#!/bin/bash
# ============================================================================
# Kspace Profile 实验封装脚本（单例 / batch 模式）
#
# 在 SimGrid 2D mesh 虚拟平台上重放真实 LAMMPS，采集应用层真实耗时的
# Kspace 操作级拆分（FFT 计算 / FFT 转置通信 / 邻居通信）。
# 方法详见 Kspace_profile.md。
#
# 用法:
#   单例模式: ./run_kspace_profile.sh -i <in.lammps> -d <data 文件> -n <rank 数>
#                                     [-s 步数] [-t 超时秒]
#   batch 模式: ./run_kspace_profile.sh -b <batch 文件> [-s 步数] [-t 超时秒]
#
# batch 文件格式（每行一个算例，# 开头为注释）:
#   <in.lammps 路径>,<data 文件路径>,[r1,r2,...][,步数]
#   例: cases/lialocl_coul_2688/in.lammps,cases/lialocl_coul_2688/data.LiAlOCl_nvt_charge,[4,8,16,32]
#        cases/lialocl_coul_1344/in.lammps,cases/lialocl_coul_1344/data.LiAlOCl_nvt_charge,[4,8],20
#   行内步数可选，缺省用全局 -s（默认 10）
#
# 单例产物 runs/<时间戳>_<np>ranks/:
#   phase_<np>.csv / run_<np>.log / log.lammps / mesh.txt   过程文件
#   result/{log.lammps,kspace_breakdown.txt}     最终结果
#
# batch 产物 runs/<时间戳>_batch/:
#   <算例名>/<np>ranks/result/{log.lammps,kspace_breakdown.txt}  各 rank 结果
#   batch_result.csv    汇总：各 rank 的 MPI breakdown + Kspace breakdown
# ============================================================================

ROOT="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
LMP=$ROOT/lammps-build-smpi/lmp
STRIP_PIE=$ROOT/simgrid_traces/strip_pie_flag.py
GEN_PLATFORM=$ROOT/simgrid_traces/gen_mesh_platform.py
PLATFORM_DIR=$ROOT/simgrid_traces
ANALYZER_DIR=$ROOT/runs/smpi_fft_split
ANALYZER=$ANALYZER_DIR/analyze_fft_split.py

STEPS=10
TIMEOUT=1500
IN_FILE=""
DATA_FILE=""
NP=""
BATCH_FILE=""

usage() {
    sed -n '2,31p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

while getopts "i:d:n:b:s:t:h" opt; do
    case $opt in
        i) IN_FILE="$OPTARG" ;;
        d) DATA_FILE="$OPTARG" ;;
        n) NP="$OPTARG" ;;
        b) BATCH_FILE="$OPTARG" ;;
        s) STEPS="$OPTARG" ;;
        t) TIMEOUT="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ "$STEPS" =~ ^[1-9][0-9]*$ ]] || { echo "错误: 步数必须是正整数: $STEPS"; exit 1; }
[ -x "$LMP" ] || { echo "错误: SMPI 版 lmp 不存在: $LMP（先按手册 4.2 编译）"; exit 1; }
command -v smpirun >/dev/null || { echo "错误: smpirun 不在 PATH"; exit 1; }

# ---------- PIE 标志检查（SMPI dlopen 要求，幂等） ----------
python3 "$STRIP_PIE" "$LMP" || { echo "错误: strip_pie_flag.py 执行失败"; exit 1; }

# ============================================================================
# 单次运行: run_one <in 文件> <data 文件> <np> <运行目录> [步数]
# 完成 平台选择/生成 -> 输入拷贝 -> smpirun -> 分析归档；失败返回非 0
# ============================================================================
run_one() {
    local in_file=$1 data_file=$2 np=$3 run_dir=$4 steps=${5:-$STEPS}

    [[ "$np" =~ ^[1-9][0-9]*$ ]] || { echo "错误: rank 数必须是正整数: $np"; return 1; }
    [[ "$steps" =~ ^[1-9][0-9]*$ ]] || { echo "错误: 运行步数必须是正整数: $steps"; return 1; }
    [ -f "$in_file" ]  || { echo "错误: in 文件不存在: $in_file"; return 1; }
    [ -f "$data_file" ] || { echo "错误: data 文件不存在: $data_file"; return 1; }

    # ---- mesh 尺寸选择：优先整除分解（KX<=KY），否则向上取整 ----
    local kx=1 ky=$np d
    for ((d = np; d >= 1; d--)); do
        if (( np % d == 0 && d * d <= np )); then
            kx=$d; ky=$((np / d)); break
        fi
    done
    if (( kx * ky != np )); then
        kx=$(python3 -c "import math; print(math.ceil(math.sqrt($np)))")
        ky=$(python3 -c "import math; print(math.ceil($np / $kx))")
    fi
    echo "      mesh 尺寸: ${kx}x${ky} ($((kx * ky)) hosts)"

    mkdir -p "$run_dir/result"
    echo "${kx}x${ky}" > "$run_dir/mesh.txt"   # 供分析/汇总读取

    # ---- 平台 XML：现成则复用，否则在运行目录生成 ----
    local platform=$PLATFORM_DIR/platform_mesh_${kx}x${ky}.xml
    if [ ! -f "$platform" ]; then
        platform=$run_dir/platform_mesh_${kx}x${ky}.xml
        python3 "$GEN_PLATFORM" "$kx" "$ky" "$run_dir" || { echo "错误: 平台生成失败"; return 1; }
    fi
    echo "      平台文件: $platform"

    # ---- 拷贝输入文件（data 按 in.lammps 中 read_data 的目标名放置） ----
    cp "$in_file" "$run_dir/in.lammps"
    local data_target
    data_target=$(grep -Em1 '^\s*read_data' "$in_file" | awk '{print $2}')
    if [ -z "$data_target" ] || [[ "$data_target" == *'$'* ]]; then
        data_target=$(basename "$data_file")
        echo "      警告: read_data 目标名无法解析，data 文件按其文件名放置: $data_target"
    fi
    mkdir -p "$run_dir/$(dirname "$data_target")"
    cp "$data_file" "$run_dir/$data_target"

    # ---- 执行 SMPI 重放 ----
    echo "      运行 smpirun: np=$np, STEPS=$steps, timeout=${TIMEOUT}s"
    (cd "$run_dir" && \
     SMPI_PRIVATIZATION=0 \
     LAMMPS_PHASE_TRACE="$run_dir/phase_${np}.csv" \
     timeout "$TIMEOUT" smpirun \
        -platform "$platform" \
        -np "$np" "$LMP" \
        -in in.lammps -v STEPS "$steps" \
        > "$run_dir/run_${np}.log" 2>&1)
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "错误: smpirun 退出码 $rc，日志: $run_dir/run_${np}.log（末尾如下）"
        tail -20 "$run_dir/run_${np}.log"
        return $rc
    fi
    local n_comm
    n_comm=$(grep -c KSPACE_FFT_COMM "$run_dir/phase_${np}.csv" 2>/dev/null || echo 0)
    echo "      完成 (rc=0), KSPACE_FFT_COMM marker: $n_comm 条"
    [ "$n_comm" -eq 0 ] && echo "      警告: 无 marker，检查输入脚本 timer 级别（需 timer normal 以上）"

    # ---- 分析：Kspace 操作拆分 ----
    python3 "$ANALYZER" "$run_dir" "$np" > "$run_dir/result/kspace_breakdown.txt" 2>&1 \
        || { echo "错误: 分析失败，输出:"; cat "$run_dir/result/kspace_breakdown.txt"; return 1; }

    # ---- 归档日志 ----
    [ -f "$run_dir/log.lammps" ] && cp "$run_dir/log.lammps" "$run_dir/result/log.lammps" \
        || echo "      警告: 未找到 log.lammps（LAMMPS 可能未正常写日志）"
    return 0
}

# ============================================================================
# batch 汇总: 遍历 <batch 目录>/<算例>/<np>ranks/，汇总 MPI breakdown
# 与 Kspace breakdown 到 batch_result.csv
# ============================================================================
summarize_batch() {
    local batch_dir=$1
    ANALYZER_DIR="$ANALYZER_DIR" python3 - "$batch_dir" <<'PYEOF'
import csv
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.environ['ANALYZER_DIR'])
from analyze_fft_split import parse_log, parse_phase, parse_mpi_breakdown, MPI_PHASES

batch = Path(sys.argv[1])
PH = [p.lower() for p in MPI_PHASES]
KS = ['fft_comm', 'fft_calc', 'neighbor', 'other']

rows = []
case_dirs = sorted(p for p in batch.iterdir() if p.is_dir())
for case_dir in case_dirs:
    run_dirs = [p for p in case_dir.iterdir()
                if p.is_dir() and re.fullmatch(r'\d+ranks', p.name)]
    run_dirs.sort(key=lambda p: int(re.match(r'\d+', p.name).group()))
    for run_dir in run_dirs:
        np = int(re.match(r'\d+', run_dir.name).group())
        row = {'case': case_dir.name, 'ranks': np}
        mesh_file = run_dir / 'mesh.txt'
        if mesh_file.exists():
            row['mesh'] = mesh_file.read_text().strip()
        log = run_dir / 'result' / 'log.lammps'
        phase_csv = run_dir / f'phase_{np}.csv'
        if not log.exists() or not phase_csv.exists():
            row['status'] = 'failed'
            rows.append(row)
            continue
        row['status'] = 'ok'
        try:
            steps, loop, _ = parse_log(log)
            row['run_steps'] = steps
            row['loop_time_s'] = f'{loop:.6f}'
            bd = parse_mpi_breakdown(log)
            for p in MPI_PHASES:
                avg, pct = bd.get(p, (None, None))
                row[f'mpi_{p.lower()}_avg_s'] = '' if avg is None else f'{avg:.6g}'
                row[f'mpi_{p.lower()}_pct_total'] = '' if pct is None else f'{pct:.2f}'
            buckets, _, ktot, _ = parse_phase(phase_csv, steps)
            for k in KS:
                row[f'ks_{k}_ms_per_step'] = f'{buckets[k] / steps * 1e3:.4f}'
                row[f'ks_{k}_pct_of_kspace'] = \
                    f'{100 * buckets[k] / ktot:.1f}' if ktot else ''
            row['ks_total_ms_per_step'] = f'{ktot / steps * 1e3:.4f}'
        except Exception as e:  # 解析异常也保留一行便于排查
            row['status'] = f'parse_error: {e}'
        rows.append(row)

fieldnames = ['case', 'ranks', 'mesh', 'status', 'run_steps', 'loop_time_s']
fieldnames += [f'mpi_{p}_avg_s' for p in PH] + [f'mpi_{p}_pct_total' for p in PH]
fieldnames += [f'ks_{k}_ms_per_step' for k in KS] + ['ks_total_ms_per_step']
fieldnames += [f'ks_{k}_pct_of_kspace' for k in KS]

out_csv = batch / 'batch_result.csv'
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, restval='')
    w.writeheader()
    w.writerows(rows)
print(f"wrote {out_csv} ({len(rows)} rows)")

# 终端摘要
print()
print(f"{'case':<24} {'ranks':>5} {'mesh':>6} {'status':>8} {'Loop(s)':>9} {'Kspace%':>8} "
      f"{'转置(ms/step)':>13} {'转置占比%':>9} {'FFT计算':>9} {'邻居通信':>9}")
for r in rows:
    print(f"{r['case']:<24} {r['ranks']:>5} {r.get('mesh', '-'):>6} {r['status']:>8} "
          f"{r.get('loop_time_s', '-'):>9} "
          f"{r.get('mpi_kspace_pct_total', '-'):>8} "
          f"{r.get('ks_fft_comm_ms_per_step', '-'):>13} "
          f"{r.get('ks_fft_comm_pct_of_kspace', '-'):>9} "
          f"{r.get('ks_fft_calc_ms_per_step', '-'):>9} "
          f"{r.get('ks_neighbor_ms_per_step', '-'):>9}")
PYEOF
}

# ============================================================================
# 模式分发
# ============================================================================
if [ -n "$BATCH_FILE" ]; then
    # ---------- batch 模式 ----------
    [ -f "$BATCH_FILE" ] || { echo "错误: batch 文件不存在: $BATCH_FILE"; exit 1; }
    BATCH_DIR=$ROOT/runs/$(date +%Y%m%d_%H%M%S)_batch
    mkdir -p "$BATCH_DIR"
    echo "batch 输出目录: $BATCH_DIR"

    N_OK=0; N_FAIL=0
    while IFS= read -r line || [ -n "$line" ]; do
        line=$(echo "$line" | sed 's/#.*//' | xargs)   # 去注释与首尾空白
        [ -z "$line" ] && continue
        in_f=$(echo "$line" | cut -d, -f1 | xargs)
        data_f=$(echo "$line" | cut -d, -f2 | xargs)
        ranks_str=$(echo "$line" | grep -o '\[[^]]*\]' | head -1 | tr -d '[] ' \
                    | sed 's/,,*/,/g; s/^,//; s/,$//')
        # 方括号后的可选步数字段，缺省用全局 -s
        line_steps=$(echo "$line" | sed -n 's/.*\][,[:space:]]*\([0-9][0-9]*\)[[:space:]]*$/\1/p')
        if [ -z "$in_f" ] || [ -z "$data_f" ] || [ -z "$ranks_str" ]; then
            echo "警告: batch 行格式错误（需 in,data,[r1,r2,...][,步数]），跳过: $line"
            N_FAIL=$((N_FAIL + 1)); continue
        fi
        # 算例子目录名取 in 文件所在目录名，重名时追加序号
        case_tag=$(basename "$(dirname "$(readlink -f "$in_f")")")
        [ "$case_tag" = "." ] || [ "$case_tag" = "/" ] && case_tag=$(basename "$in_f" .lammps)
        tag=$case_tag; i=2
        while [ -d "$BATCH_DIR/$tag" ]; do tag=${case_tag}_$i; i=$((i + 1)); done

        echo "=============================================================="
        echo "算例: $tag  ($in_f)  步数: ${line_steps:-$STEPS}"
        IFS=',' read -ra RANK_LIST <<< "$ranks_str"
        for np in "${RANK_LIST[@]}"; do
            echo "--------------------------------------------------------------"
            echo ">>> ranks=$np"
            if run_one "$in_f" "$data_f" "$np" "$BATCH_DIR/$tag/${np}ranks" "$line_steps"; then
                N_OK=$((N_OK + 1))
            else
                echo ">>> ranks=$np 失败，继续后续配置"
                N_FAIL=$((N_FAIL + 1))
            fi
        done
    done < "$BATCH_FILE"

    echo "=============================================================="
    echo "全部运行结束: 成功 $N_OK, 失败 $N_FAIL。生成汇总 CSV..."
    summarize_batch "$BATCH_DIR" || { echo "错误: 汇总失败"; exit 1; }
    echo
    echo "batch 结果: $BATCH_DIR/batch_result.csv"
    [ $N_FAIL -gt 0 ] && exit 2
    exit 0
else
    # ---------- 单例模式 ----------
    [ -z "$IN_FILE" ] || [ -z "$DATA_FILE" ] || [ -z "$NP" ] && usage
    RUN_DIR=$ROOT/runs/$(date +%Y%m%d_%H%M%S)_${NP}ranks
    echo "运行目录: $RUN_DIR"
    run_one "$IN_FILE" "$DATA_FILE" "$NP" "$RUN_DIR" || exit $?
    echo "完成。结果:"
    echo "      $RUN_DIR/result/log.lammps"
    echo "      $RUN_DIR/result/kspace_breakdown.txt"
    echo
    grep -A20 "趋势表（ms/step，真实时钟）" "$RUN_DIR/result/kspace_breakdown.txt" || \
        cat "$RUN_DIR/result/kspace_breakdown.txt"
fi
