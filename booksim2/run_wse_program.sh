#!/usr/bin/env bash
# Run a Phase-1 WSE program through BookSim's manager-level multicast replay.
set -euo pipefail

BS_DIR="$(cd "$(dirname "$0")" && pwd)"
BOOKSIM="$BS_DIR/src/booksim"
TEMPLATE="${2:-$BS_DIR/ccdg_lammps_4x4.cfg}"
PROGRAM="${1:?usage: $0 <wse_phase1.program.json> [booksim-template.cfg] [timeout]}"
TIMEOUT="${3:-3600}"
PROGRAM="$(cd "$(dirname "$PROGRAM")" && pwd)/$(basename "$PROGRAM")"
PREFIX="${PROGRAM%.program.json}"
REPLAY="$PREFIX.replay.ccdg"
EST="$PREFIX.replay.est"
REPORT="$PREFIX.report.json"
CFG="$PREFIX.booksim.cfg"
LOG="$PREFIX.booksim.log"
STATS="$PREFIX.booksim.stats"
ACCEPTANCE="$PREFIX.acceptance.json"

for file in "$BOOKSIM" "$TEMPLATE" "$PROGRAM" "$REPLAY" "$EST" "$REPORT"; do
  [ -e "$file" ] || { echo "ERROR: missing WSE replay input $file" >&2; exit 2; }
done

read -r K FLIT CAP EXPECTED PROFILE <<<"$(python3 - "$PROGRAM" "$REPORT" <<'PY'
import json, sys
p, r = (json.load(open(path)) for path in sys.argv[1:])
mesh = p["hardware"]["mesh"]
if mesh[0] != mesh[1]:
    raise SystemExit("WSE replay requires a square mesh")
print(mesh[0], p["hardware"]["flit_size_bytes"],
      p["hardware"]["compute_capability_ops_s"], r["cycles"]["compiled_total"],
      p["hardware"].get("profile", "cfg_native"))
PY
)"

if [ "$PROFILE" = wse_fast ]; then
  PROFILE_SED=(-e "s|^num_vcs = .*|num_vcs = 4;|"
               -e "s|^routing_delay = .*|routing_delay = 0;|"
               -e "s|^vc_alloc_delay = .*|vc_alloc_delay = 1;|"
               -e "s|^sw_alloc_delay = .*|sw_alloc_delay = 1;|")
else
  PROFILE_SED=()
fi
sed -e "s|^sim_type = .*|sim_type = wse;|" \
    -e "s|^ccdg_file = .*|ccdg_file = $REPLAY;|" \
    -e "s|^ccdg_schedule_file = .*|ccdg_schedule_file = $EST;|" \
    -e "s|^flit_size_bytes = .*|flit_size_bytes = $FLIT;|" \
    -e "s|^ccdg_compute_capability = .*|ccdg_compute_capability = $CAP;|" \
    -e "s|^k = .*|k = $K;|" \
    -e "s|^wse_gating = .*|wse_gating = 0;|" \
    -e "s|^stats_out = .*|stats_out = $STATS;|" \
    "${PROFILE_SED[@]}" \
    "$TEMPLATE" > "$CFG"
printf '\nwse_program_file = %s;\n' "$PROGRAM" >> "$CFG"

rc=0
timeout "$TIMEOUT" "$BOOKSIM" "$CFG" > "$LOG" 2>&1 || rc=$?
[ "$rc" -eq 0 ] || [ "$rc" -eq 255 ] || {
  echo "ERROR: BookSim WSE replay failed rc=$rc; see $LOG" >&2
  exit "$rc"
}
MEASURED="$(grep -oP 'completed in \K[0-9]+' "$LOG" | tail -1 || true)"
[ -n "$MEASURED" ] || { echo "ERROR: no BookSim completion cycle in $LOG" >&2; exit 1; }

python3 - "$EXPECTED" "$MEASURED" "$PROGRAM" "$CFG" "$LOG" "$STATS" "$ACCEPTANCE" <<'PY'
import json, re, sys
expected, measured = map(int, sys.argv[1:3])
error = abs(measured - expected) / expected if expected else float("inf")
stats_text = open(sys.argv[6]).read()
def stat(name, cast=int):
    match = re.search(rf"^{name} = ([0-9.eE+-]+);", stats_text, re.M)
    if not match:
        raise SystemExit(f"missing WSE stat: {name}")
    return cast(match.group(1))
wave_expected = stat("wse_wavefronts_expected")
wave_injected = stat("wse_wavefronts_injected")
branch_expected = stat("wse_branches_expected")
branch_delivered = stat("wse_branches_delivered")
command_expected = stat("wse_commands_expected")
command_delivered = stat("wse_commands_delivered")
congestion = stat("wse_congestion_ratio", float)
average_packet_queue_cycles = stat("average_packet_queue_cycles", float)
average_flit_queue_cycles = stat("average_flit_queue_cycles", float)
average_injection_rate = stat("average_injection_rate", float)
injection_saturation_ratio = stat("injection_saturation_ratio", float)
saturated_injection_rate = stat("saturated_injection_rate", float)
communication_to_compute_ratio = stat("communication_to_compute_ratio", float)
exposed_communication_to_compute_ratio = stat("exposed_communication_to_compute_ratio", float)
checks = {
    "cycles_within_1pct": error <= 0.01,
    "wavefronts_conserved": wave_injected == wave_expected,
    "branches_conserved": branch_delivered == branch_expected,
    "commands_conserved": command_delivered == command_expected,
    "congestion_le_0_1pct": congestion <= 0.001,
}
result = {
    "status": "PASS" if all(checks.values()) else "FAIL",
    "compiler_cycles": expected,
    "booksim_cycles": measured,
    "relative_error": error,
    "threshold": 0.01,
    "wavefronts": {"expected": wave_expected, "injected": wave_injected},
    "branches": {"expected": branch_expected, "delivered": branch_delivered},
    "commands": {"expected": command_expected, "delivered": command_delivered},
    "congestion_ratio": congestion,
    "average_packet_queue_cycles": average_packet_queue_cycles,
    "average_flit_queue_cycles": average_flit_queue_cycles,
    "average_injection_rate": average_injection_rate,
    "injection_saturation_ratio": injection_saturation_ratio,
    "saturated_injection_rate": saturated_injection_rate,
    "communication_to_compute_ratio": communication_to_compute_ratio,
    "exposed_communication_to_compute_ratio": exposed_communication_to_compute_ratio,
    "checks": checks,
    "program": sys.argv[3],
    "booksim_cfg": sys.argv[4],
    "booksim_log": sys.argv[5],
}
with open(sys.argv[7], "w") as fh:
    json.dump(result, fh, indent=2)
    fh.write("\n")
print("WSE BookSim acceptance: "
      f"{result['status']} compiler={expected} measured={measured} error={error:.4%} "
      f"wavefronts={wave_injected}/{wave_expected} "
      f"branches={branch_delivered}/{branch_expected} congestion={congestion:.4%} "
      f"queue={average_packet_queue_cycles:.3f}cyc "
      f"sat_injection={saturated_injection_rate:.4%} "
      f"comm/compute={communication_to_compute_ratio:.6g}")
raise SystemExit(0 if result["status"] == "PASS" else 1)
PY
