#!/usr/bin/env bash
# Per-rank hardware DP operation counter for LAMMPS.
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
    -e "$FLOPS_EVENTS" -- "$@"
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
FLOPS_OUT_DIR=${FLOPS_OUT_DIR:-"$(dirname "$INPUT")/flops_profile"}
FLOPS_REPEATS=${FLOPS_REPEATS:-3}

[[ "$RANKS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid ranks: $RANKS" >&2; exit 2; }
[[ "$STEPS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid steps: $STEPS" >&2; exit 2; }
[[ "$FLOPS_REPEATS" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: invalid FLOPS_REPEATS: $FLOPS_REPEATS" >&2
  exit 2
}
[ -f "$INPUT" ] || { echo "ERROR: missing input: $INPUT" >&2; exit 2; }
[ -x "$LMP" ] || { echo "ERROR: missing LAMMPS binary: $LMP" >&2; exit 2; }
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
if probe_events "${NAMED_FP_EVENTS[@]}"; then
  FP_EVENTS=("${NAMED_FP_EVENTS[@]}")
  EVENT_ENCODING=named
elif [ "$(uname -m)" = x86_64 ] && probe_events "${RAW_FP_EVENTS[@]}"; then
  FP_EVENTS=("${RAW_FP_EVENTS[@]}")
  EVENT_ENCODING=intel_raw_c7
else
  echo "ERROR: Intel FP_ARITH_INST_RETIRED counters are unavailable" >&2
  cat "$probe" >&2
  exit 2
fi
EVENTS=(cycles instructions "${FP_EVENTS[@]}")
FLOPS_EVENTS=$(IFS=,; echo "${EVENTS[*]}")

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
  export PERF_BIN FLOPS_OUT_DIR FLOPS_EVENTS FLOPS_TAG=$1
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

python3 - "$FLOPS_OUT_DIR" "$RANKS" "$STEPS" "$FLOPS_REPEATS" \
  "$EVENT_ENCODING" "${EVENTS[@]}" <<'PY'
import csv
import json
import pathlib
import statistics
import sys

out = pathlib.Path(sys.argv[1])
ranks, steps = map(int, sys.argv[2:4])
repeats = int(sys.argv[4])
encoding = sys.argv[5]
events = sys.argv[6:]
weights = {
    "fp_arith_inst_retired.scalar_double": 1,
    "fp_arith_inst_retired.128b_packed_double": 2,
    "fp_arith_inst_retired.256b_packed_double": 4,
    "fp_arith_inst_retired.512b_packed_double": 8,
    "r01c7": 1,
    "r04c7": 2,
    "r10c7": 4,
    "r40c7": 8,
}

def read(path):
    values = {}
    with path.open() as f:
        for row in csv.reader(f, delimiter=";"):
            event = row[2].strip().split(":", 1)[0] if len(row) >= 3 else ""
            if event not in events:
                continue
            value = row[0].strip().replace(",", "")
            if value.startswith("<"):
                raise RuntimeError(f"{path}: {event} = {value}")
            values[event] = float(value)
    missing = set(events) - values.keys()
    if missing:
        raise RuntimeError(f"{path}: missing events: {sorted(missing)}")
    return values

records = []
for rank in range(ranks):
    samples = []
    for repeat in range(repeats):
        base = read(out / f"baseline.repeat{repeat}.rank{rank:04d}.csv")
        run = read(out / f"run.repeat{repeat}.rank{rank:04d}.csv")
        delta = {event: run[event] - base[event] for event in events}
        dp_ops = sum(delta[event] * weights[event] for event in events[2:])
        samples.append({
            "repeat": repeat,
            "dp_ops_per_step": dp_ops / steps,
            "cycles_per_step": delta["cycles"] / steps,
            "instructions_per_step": delta["instructions"] / steps,
            "baseline_events": base,
            "run_events": run,
            "delta_events": delta,
        })
    records.append({
        "rank": rank,
        "steps": steps,
        "repeats": repeats,
        "aggregation": "per-rank median",
        "dp_ops_per_step": statistics.median(s["dp_ops_per_step"] for s in samples),
        "cycles_per_step": statistics.median(s["cycles_per_step"] for s in samples),
        "instructions_per_step": statistics.median(
            s["instructions_per_step"] for s in samples
        ),
        "repeat_records": samples,
    })

per_step = [r["dp_ops_per_step"] for r in records]
summary = {
    "schema_version": 1,
    "method": (
        "median of repeated perf [run(N)-run(0)]/N, "
        "Intel FP_ARITH_INST_RETIRED lane weighted"
    ),
    "event_encoding": encoding,
    "event_semantics": (
        "Intel event 0xC7 double-precision arithmetic, weighted by SIMD lanes; "
        "raw umasks: 0x01 scalar, 0x04 128b, 0x10 256b, 0x40 512b"
    ),
    "num_ranks": ranks,
    "steps": steps,
    "repeats": repeats,
    "dp_ops_per_step": {
        "sum": sum(per_step),
        "average": sum(per_step) / ranks,
        "minimum": min(per_step),
        "maximum": max(per_step),
    },
    "ranks": records,
}
(out / "flops_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary["dp_ops_per_step"], sort_keys=True))
PY

echo "FLOPS profile: $FLOPS_OUT_DIR/flops_summary.json"
