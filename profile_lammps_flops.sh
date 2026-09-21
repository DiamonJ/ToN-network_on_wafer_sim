#!/usr/bin/env bash
# Per-rank instruction, cycle, and wall-time profiler for LAMMPS.
# Usage: ./profile_lammps_flops.sh RANKS in.lammps [STEPS=100] [LMP=install/bin/lmp]
set -euo pipefail

ROOT="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
SELF=$(readlink -f "${BASH_SOURCE[0]}")

if [ "${1:-}" = --rank ]; then
  shift
  rank=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-${PMIX_RANK:-}}}
  [ -n "$rank" ] || { echo "ERROR: cannot determine MPI rank" >&2; exit 2; }
  printf -v rank "%04d" "$rank"
  exec "$PERF_BIN" stat -x ';' -o "$FLOPS_OUT_DIR/$FLOPS_TAG.rank${rank}.csv" \
    -e "$PROFILE_EVENTS" -- "$@"
fi

[ "$#" -ge 2 ] || {
  echo "Usage: $0 RANKS in.lammps [STEPS=100] [LMP=install/bin/lmp]" >&2
  exit 2
}

RANKS=$1
INPUT=$(readlink -f "$2")
STEPS=${3:-100}
LMP=$(readlink -f "${4:-$ROOT/install/bin/lmp}")
PERF_BIN=${PERF_BIN:-perf}
FLOPS_OUT_DIR=${COMPUTE_PROFILE_DIR:-${FLOPS_OUT_DIR:-"$(dirname "$INPUT")/compute_profile"}}
FLOPS_REPEATS=${COMPUTE_REPEATS:-${FLOPS_REPEATS:-3}}
SUMMARIZER=$ROOT/summarize_compute_profile.py

[[ "$RANKS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid ranks: $RANKS" >&2; exit 2; }
[[ "$STEPS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid steps: $STEPS" >&2; exit 2; }
[[ "$FLOPS_REPEATS" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: invalid FLOPS_REPEATS: $FLOPS_REPEATS" >&2
  exit 2
}
[ -f "$INPUT" ] || { echo "ERROR: missing input: $INPUT" >&2; exit 2; }
[ -x "$LMP" ] || { echo "ERROR: missing LAMMPS binary: $LMP" >&2; exit 2; }
[ -f "$SUMMARIZER" ] || { echo "ERROR: missing profile summarizer: $SUMMARIZER" >&2; exit 2; }
command -v mpirun >/dev/null || { echo "ERROR: mpirun not found" >&2; exit 2; }
command -v "$PERF_BIN" >/dev/null || {
  echo "ERROR: perf not found; install the linux-tools/perf package for this kernel" >&2
  exit 2
}

NAMED_FP_EVENTS=(
  fp_arith_inst_retired.scalar_double
  fp_arith_inst_retired.128b_packed_double
  fp_arith_inst_retired.256b_packed_double
  fp_arith_inst_retired.512b_packed_double
)
# Intel event 0xC7: umask 0x01 is scalar double (0x02 is scalar single).
RAW_FP_EVENTS=(r01c7 r04c7 r10c7 r40c7)

# Prefer readable aliases; old/custom kernels may require Intel raw encodings.
probe=$(mktemp)
trap 'rm -f "$probe"' EXIT
probe_events() {
  local event
  for event in "$@"; do
    "$PERF_BIN" stat -x ';' -o "$probe" -e "$event" -- \
      python3 -c 'sum(i * i for i in range(100000))' \
      2>/dev/null || return 1
    grep -Eq '<not supported>|<not counted>' "$probe" && return 1
  done
  return 0
}
EVENTS=(cycles instructions duration_time)
EVENT_ENCODING=core
if probe_events "${NAMED_FP_EVENTS[@]}"; then
  FP_EVENTS=("${NAMED_FP_EVENTS[@]}")
  EVENT_ENCODING=core+named_fp
elif [ "$(uname -m)" = x86_64 ] && probe_events "${RAW_FP_EVENTS[@]}"; then
  FP_EVENTS=("${RAW_FP_EVENTS[@]}")
  EVENT_ENCODING=core+intel_raw_c7
else
  FP_EVENTS=()
  echo "WARN: FP_ARITH_INST_RETIRED unavailable; continuing with instructions/cycles/time" >&2
fi
EVENTS+=("${FP_EVENTS[@]}")
PROFILE_EVENTS=$(IFS=,; echo "${EVENTS[*]}")

read_ghz() { awk -v khz="$(cat "$1")" 'BEGIN { printf "%.9g", khz / 1000000.0 }'; }
if [ -n "${CPU_FREQUENCY_GHZ:-}" ]; then
  CPU_FREQ_NOMINAL_GHZ=$CPU_FREQUENCY_GHZ
  CPU_FREQ_MIN_GHZ=${CPU_FREQUENCY_MIN_GHZ:-$CPU_FREQUENCY_GHZ}
  CPU_FREQ_MAX_GHZ=${CPU_FREQUENCY_MAX_GHZ:-$CPU_FREQUENCY_GHZ}
  CPU_FREQ_SOURCE=environment
elif [ -r /sys/devices/system/cpu/cpu0/cpufreq/base_frequency ]; then
  CPU_FREQ_NOMINAL_GHZ=$(read_ghz /sys/devices/system/cpu/cpu0/cpufreq/base_frequency)
  CPU_FREQ_MIN_GHZ=$(read_ghz /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_min_freq)
  CPU_FREQ_MAX_GHZ=$(read_ghz /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq)
  CPU_FREQ_SOURCE=sysfs_cpufreq
else
  read -r CPU_FREQ_MIN_GHZ CPU_FREQ_NOMINAL_GHZ CPU_FREQ_MAX_GHZ < <(
    awk -F: '/cpu MHz/{gsub(/ /,"",$2); a[++n]=$2/1000; s+=$2/1000}
      END {if (!n) exit 1; min=a[1]; max=a[1]; for(i=2;i<=n;i++){if(a[i]<min)min=a[i];if(a[i]>max)max=a[i]} print min,s/n,max}' \
      /proc/cpuinfo
  )
  CPU_FREQ_SOURCE=proc_cpuinfo_snapshot
fi
awk -v lo="$CPU_FREQ_MIN_GHZ" -v mid="$CPU_FREQ_NOMINAL_GHZ" -v hi="$CPU_FREQ_MAX_GHZ" \
  'BEGIN { exit !(lo > 0 && lo <= mid && mid <= hi) }' || {
    echo "ERROR: CPU frequency must satisfy 0 < min <= nominal <= max" >&2
    exit 2
  }

mkdir -p "$FLOPS_OUT_DIR"
WORKDIR=$(dirname "$INPUT")
BASE_INPUT=$FLOPS_OUT_DIR/in.run0.lammps
RUN_INPUT=$FLOPS_OUT_DIR/in.run${STEPS}.lammps

replace_last_run() {
  awk -v steps="$2" '
    { line[NR]=$0; if ($1 == "run") last=NR }
    END {
      if (!last) exit 2
      sub(/^[[:space:]]*run[[:space:]].*/, "run             " steps, line[last])
      for (i=1; i<=NR; i++) print line[i]
    }' "$1"
}
replace_last_run "$INPUT" 0 > "$BASE_INPUT"
replace_last_run "$INPUT" "$STEPS" > "$RUN_INPUT"

run_profile() {
  export PERF_BIN FLOPS_OUT_DIR PROFILE_EVENTS FLOPS_TAG=$1
  (
    cd "$WORKDIR"
    unset LD_PRELOAD DUMPI_OUTDIR LAMMPS_WSE_PLAN
    mpirun -np "$RANKS" --allow-run-as-root --oversubscribe \
      "$SELF" --rank "$LMP" -screen none -log none -in "$2"
  )
}

for ((repeat=0; repeat<FLOPS_REPEATS; repeat++)); do
  run_profile "baseline.repeat${repeat}" "$BASE_INPUT"
  run_profile "run.repeat${repeat}" "$RUN_INPUT"
done

python3 "$SUMMARIZER" "$FLOPS_OUT_DIR" "$RANKS" "$STEPS" "$FLOPS_REPEATS" \
  "$CPU_FREQ_NOMINAL_GHZ" "$CPU_FREQ_MIN_GHZ" "$CPU_FREQ_MAX_GHZ" \
  "$CPU_FREQ_SOURCE" "$EVENT_ENCODING" "${EVENTS[@]}" \
  --input-file "$INPUT" \
  -o "$FLOPS_OUT_DIR/compute_profile.json"

echo "Compute profile: $FLOPS_OUT_DIR/compute_profile.json"
