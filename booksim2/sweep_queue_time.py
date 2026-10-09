#!/usr/bin/env python3
"""
BookSim injection-rate sweep → average queueing time plot.

Runs standard BookSim with synthetic traffic patterns,
sweeps injection_rate, and plots:
  X = injection_rate (flits/node/cycle)
  Y = average queueing time (avg_plat - avg_nlat, in cycles)

Uses "include_queuing=1" so that queueing time = head_flit waiting
time in the source injection queue before entering the network.

Scenarios:
  1. 4x4 mesh, uniform traffic
  2. 8x8 mesh, uniform traffic
  3. 4x4 mesh, transpose traffic
  4. 4x4 mesh, uniform traffic + WSE gating
"""

import subprocess
import sys
import os
import csv
import json
from pathlib import Path
from typing import Optional, List, Tuple

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BOOKSIM = Path(__file__).resolve().parent / "booksim"
if not BOOKSIM.exists():
    BOOKSIM = Path(__file__).resolve().parent / "src" / "booksim"

OUTPUT_DIR = Path(__file__).resolve().parent / "sweep_results"

# Injection rate range
INJ_START = 0.01
INJ_STEP  = 0.005
INJ_STOP  = 0.61   # exclusive stop; sweeps normally start failing around 0.3-0.5

# Per-run settings
SIM_COUNT = 1
WARMUP_PERIODS = 3
SAMPLE_PERIOD  = 1000
MAX_SAMPLES    = 50
LATENCY_THRES  = 10000.0
PACKET_SIZE    = 8
NUM_VCS        = 24
VC_BUF_SIZE    = 8
INCLUDE_QUEUING = 1     # so plat - nlat gives injection queue delay

# ---------------------------------------------------------------------------
# Scenario definitions
# ---------------------------------------------------------------------------
SCENARIOS = [
    {
        "label": "4×4 mesh, uniform",
        "extra": "k=4 n=2 traffic=uniform wse_gating=0",
    },
    {
        "label": "8×8 mesh, uniform",
        "extra": "k=8 n=2 traffic=uniform wse_gating=0",
    },
    {
        "label": "4×4 mesh, transpose",
        "extra": "k=4 n=2 traffic=transpose wse_gating=0",
    },
    {
        "label": "4×4 mesh, uniform + WSE gating",
        "extra": "k=4 n=2 traffic=uniform wse_gating=1 wse_phase_width=64 wse_strip_width=2 wse_vc_split=1",
    },
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def run_booksim(injection_rate: float, extra_params: str) -> Optional[float]:
    """Run BookSim once.  Return average queueing time (cycles) or None on failure."""
    # Use the existing traffic_wse_4x4.cfg as a base config;
    # its k=4,n=2 will be overridden by extra_params.
    cfg_path = Path(__file__).resolve().parent / "traffic_wse_4x4.cfg"
    cmd = (
        f"{BOOKSIM} {cfg_path} "
        f"injection_rate={injection_rate} "
        f"print_csv_results=1 "
        f"stats_out=/dev/null "
        f"include_queuing={INCLUDE_QUEUING} "
        f"{extra_params}"
    )
    try:
        proc = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=120
        )
    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT] injection_rate={injection_rate}")
        return None

    for line in proc.stdout.splitlines():
        if line.startswith("results:"):
            parts = line.split(",")
            # Field layout (0-indexed):
            #   0=class  1=traffic  2=use_read_write  3=inj_rate
            #   4=min_plat  5=avg_plat  6=max_plat
            #   7=min_nlat  8=avg_nlat  9=max_nlat  ...
            if len(parts) < 9:
                return None
            try:
                avg_plat = float(parts[5])
                avg_nlat = float(parts[8])
            except ValueError:
                return None
            queue_time = avg_plat - avg_nlat
            return queue_time

    # If no results line found, the sim may have been unstable
    return None


def sweep_scenario(label: str, extra: str, start: float, step: float, stop: float) -> List[Tuple[float, float]]:
    """Sweep injection_rate for one scenario.  Return list of (inj, qtime)."""
    print(f"\n{'='*60}")
    print(f"Scenario: {label}")
    print(f"{'='*60}")
    points = []
    inj = start
    while inj < stop:
        qtime = run_booksim(inj, extra)
        if qtime is not None:
            print(f"  inj={inj:.4f}  queue_time={qtime:.2f}")
            points.append((inj, qtime))
        else:
            print(f"  inj={inj:.4f}  [FAILED / UNSTABLE] — stopping sweep")
            break
        inj = round(inj + step, 10)
    return points


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    all_data = {}

    for sc in SCENARIOS:
        points = sweep_scenario(
            sc["label"], sc["extra"],
            INJ_START, INJ_STEP, INJ_STOP,
        )
        all_data[sc["label"]] = points

        # Save CSV
        csv_path = OUTPUT_DIR / f"{sc['label'].replace(' ', '_').replace('×','x')}.csv"
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["injection_rate", "avg_queue_time_cycles"])
            w.writerows(points)
        print(f"  → saved {csv_path}")

    # Save all_data.json for plotting
    json_path = OUTPUT_DIR / "all_data.json"
    with open(json_path, "w") as f:
        json.dump(all_data, f, indent=2)
    print(f"\nAll data saved → {json_path}")

    # ---------- plotting ----------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not available; skipping plot. Data saved to CSV/JSON.")
        return

    plt.figure(figsize=(10, 6))

    markers = ["o", "s", "^", "D"]
    for i, (label, points) in enumerate(all_data.items()):
        if not points:
            continue
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        plt.plot(xs, ys, marker=markers[i % len(markers)], markersize=5,
                 linewidth=1.5, label=label)

    plt.xlabel("Injection Rate (flits/node/cycle)", fontsize=12)
    plt.ylabel("Average Queueing Time (cycles)", fontsize=12)
    plt.title("Injection Rate vs Average Queueing Time", fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    png_path = OUTPUT_DIR / "inj_vs_queue_time.png"
    plt.savefig(png_path, dpi=150)
    print(f"Plot saved → {png_path}")

    svg_path = OUTPUT_DIR / "inj_vs_queue_time.svg"
    plt.savefig(svg_path)
    print(f"Plot saved → {svg_path}")


if __name__ == "__main__":
    main()