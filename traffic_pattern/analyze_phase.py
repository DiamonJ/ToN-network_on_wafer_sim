#!/usr/bin/env python3
"""
analyze_phase.py — marker-driven global/local communication attribution.

Correlates a CCDG (with per-node wall_time_sec / wall_duration_sec emitted by
dumpi2ccdg) with a LAMMPS phase trace (phase_trace.csv emitted by the
instrumented Timer, CLOCK_MONOTONIC ns) to measure, WITHOUT any message-size
heuristic:

  * how much communication wall time falls in each top-level LAMMPS phase
  * inside the Kspace phase: the split between
      - LOCAL  communication (PPPM grid ghost exchange:
                KSPACE_GRID_REVERSE / KSPACE_GRID_FORWARD)
      - GLOBAL communication (FFT transpose remaps + global reductions:
                KSPACE_FFT / KSPACE_REDUCE)

phase_trace.csv semantics (see timer.cpp):
  - rows like  "<ns>,Kspace"  mark the END of a timer phase: the interval
    (prev_ts, ts] belongs to the named phase.
  - rows like  "<ns>,KSPACE_FFT" are START markers of sub-stages inside the
    Kspace phase: the segment [ts, next_marker_ts) belongs to the sub-stage.

CLOCK_MONOTONIC is system-wide, so rank-0 phase windows are valid for the
wall_time_sec of every rank in the CCDG.

Usage:
  single: python3 analyze_phase.py <ccdg.json> <phase_trace.csv> [--json OUT]
  batch : python3 analyze_phase.py --batch <run_dir> [<run_dir> ...]
                                          [--csv OUT] [--json OUT]
"""

import argparse
import bisect
import csv
import glob
import json
import os
import sys

# Kspace sub-phase -> communication class
SUBPHASE_CLASS = {
    "KSPACE_GRID_REVERSE": "kspace_local",
    "KSPACE_GRID_FORWARD": "kspace_local",
    "KSPACE_FFT":          "kspace_global",
    "KSPACE_REDUCE":       "kspace_global",
    "KSPACE_FORCE":        "kspace_compute",
}
GLOBAL_SUBPHASES = ("KSPACE_FFT", "KSPACE_REDUCE")
LOCAL_SUBPHASES = ("KSPACE_GRID_REVERSE", "KSPACE_GRID_FORWARD")


# ----------------------------------------------------------------------
# phase trace parsing
# ----------------------------------------------------------------------
def parse_phase_trace(path):
    """Return (top_intervals, kspace_windows).

    top_intervals : list of (end_ns, phase_name), one per top-level phase end
    kspace_windows: list of dicts {start_ns, end_ns, markers:[(ns,name)]}
    """
    events = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",", 1)
            if len(parts) != 2:
                continue
            try:
                ts = int(parts[0])
            except ValueError:
                continue
            events.append((ts, parts[1].strip()))
    events.sort(key=lambda e: e[0])

    top_intervals = [(ts, name) for ts, name in events
                     if not name.startswith("KSPACE_")]

    # top-level end timestamps for interval start lookup
    top_ts = [ts for ts, _ in top_intervals]

    kspace_windows = []
    for i, (ts, name) in enumerate(top_intervals):
        if name != "Kspace":
            continue
        start = top_ts[i - 1] if i > 0 else 0
        markers = [(t, n) for t, n in events
                   if start < t < ts and n.startswith("KSPACE_")]
        kspace_windows.append({"start": start, "end": ts, "markers": markers})

    return top_intervals, kspace_windows


# ----------------------------------------------------------------------
# classification
# ----------------------------------------------------------------------
class Classifier:
    def __init__(self, top_intervals, kspace_windows):
        self.top = top_intervals
        self.top_ends = [ts for ts, _ in top_intervals]
        self.win_starts = [w["start"] for w in kspace_windows]
        self.windows = kspace_windows

    def top_phase(self, t_ns):
        """Top-level phase owning time t_ns (end-marker semantics)."""
        i = bisect.bisect_left(self.top_ends, t_ns)
        if i >= len(self.top):
            return "outside"
        return self.top[i][1]

    def classify(self, t_ns):
        """Return (top_phase, kspace_subphase_or_None)."""
        wi = bisect.bisect_right(self.win_starts, t_ns) - 1
        if wi >= 0:
            w = self.windows[wi]
            if t_ns <= w["end"]:
                # inside a Kspace window: locate sub-stage (start markers)
                mk_ts = [m[0] for m in w["markers"]]
                mi = bisect.bisect_right(mk_ts, t_ns) - 1
                if mi < 0:
                    return "Kspace", "KSPACE_PRE"
                return "Kspace", w["markers"][mi][1]
        return self.top_phase(t_ns), None


# ----------------------------------------------------------------------
# CCDG analysis
# ----------------------------------------------------------------------
def analyze(ccdg_path, phase_csv_path):
    with open(ccdg_path) as f:
        ccdg = json.load(f)
    top_intervals, kspace_windows = parse_phase_trace(phase_csv_path)
    clf = Classifier(top_intervals, kspace_windows)

    num_ranks = ccdg.get("num_ranks", 0)
    t0 = phase_csv_path

    # accumulators
    phase_acc = {}          # phase -> [events, bytes, duration]
    sub_acc = {}            # kspace subphase -> [events, bytes, duration]
    total = {"events": 0, "bytes": 0, "dur": 0.0}
    missing_wall = 0

    for node in ccdg["nodes"]:
        if node["type"] == "COMPUTE":
            continue
        wt = node.get("wall_time_sec")
        if wt is None or wt < 0:
            missing_wall += 1
            continue
        dur = node.get("wall_duration_sec", 0.0)
        nbytes = node.get("comm_bytes", 0)
        t_ns = int(wt * 1e9)
        phase, sub = clf.classify(t_ns)

        a = phase_acc.setdefault(phase, [0, 0, 0.0])
        a[0] += 1; a[1] += nbytes; a[2] += dur
        if phase == "Kspace" and sub is not None:
            s = sub_acc.setdefault(sub, [0, 0, 0.0])
            s[0] += 1; s[1] += nbytes; s[2] += dur
        total["events"] += 1; total["bytes"] += nbytes; total["dur"] += dur

    kspace = phase_acc.get("Kspace", [0, 0, 0.0])
    g = [0, 0, 0.0]; l = [0, 0, 0.0]
    for name in GLOBAL_SUBPHASES:
        s = sub_acc.get(name, [0, 0, 0.0])
        g[0] += s[0]; g[1] += s[1]; g[2] += s[2]
    for name in LOCAL_SUBPHASES:
        s = sub_acc.get(name, [0, 0, 0.0])
        l[0] += s[0]; l[1] += s[1]; l[2] += s[2]

    def pct(x, y):
        return 100.0 * x / y if y > 0 else 0.0

    result = {
        "ccdg": os.path.basename(ccdg_path),
        "phase_trace": os.path.basename(t0),
        "num_ranks": num_ranks,
        "num_kspace_windows": len(kspace_windows),
        "missing_wall_time_nodes": missing_wall,
        "total_comm": {
            "events": total["events"], "bytes": total["bytes"],
            "wall_sec": total["dur"],
        },
        "phase_breakdown": {
            p: {"events": v[0], "bytes": v[1], "wall_sec": v[2],
                "pct_of_comm_wall": pct(v[2], total["dur"])}
            for p, v in sorted(phase_acc.items(), key=lambda kv: -kv[1][2])
        },
        "kspace_subphase_breakdown": {
            s: {"events": v[0], "bytes": v[1], "wall_sec": v[2],
                "pct_of_kspace_wall": pct(v[2], kspace[2])}
            for s, v in sorted(sub_acc.items(), key=lambda kv: -kv[1][2])
        },
        "summary": {
            "kspace_pct_of_comm_wall": pct(kspace[2], total["dur"]),
            "kspace_global_wall_sec": g[2],
            "kspace_local_wall_sec": l[2],
            "kspace_global_events": g[0],
            "kspace_local_events": l[0],
            "kspace_global_bytes": g[1],
            "kspace_local_bytes": l[1],
            "global_pct_of_kspace_wall": pct(g[2], kspace[2]),
            "global_pct_of_comm_wall": pct(g[2], total["dur"]),
            "local_pct_of_kspace_wall": pct(l[2], kspace[2]),
            "nonkspace_pct_of_comm_wall": pct(total["dur"] - kspace[2],
                                              total["dur"]),
        },
    }
    return result


# ----------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------
def print_report(r):
    s = r["summary"]
    print("=" * 72)
    print(f"Phase communication analysis: ranks={r['num_ranks']}  "
          f"({r['ccdg']} + {r['phase_trace']})")
    print("=" * 72)
    t = r["total_comm"]
    print(f"Total communication: {t['events']} events, "
          f"{t['bytes']/1e6:.2f} MB, wall {t['wall_sec']:.4f} s"
          f"   (kspace windows: {r['num_kspace_windows']})")
    if r["missing_wall_time_nodes"]:
        print(f"  WARNING: {r['missing_wall_time_nodes']} comm nodes lack "
              f"wall_time_sec and were excluded")
    print("\n--- Communication wall time by LAMMPS phase ---")
    for p, v in r["phase_breakdown"].items():
        print(f"  {p:<10s} {v['wall_sec']:10.4f} s  "
              f"{v['pct_of_comm_wall']:6.2f}%   "
              f"({v['events']} events, {v['bytes']/1e6:.2f} MB)")
    print("\n--- Kspace sub-phase breakdown (marker-driven) ---")
    for sp, v in r["kspace_subphase_breakdown"].items():
        cls = SUBPHASE_CLASS.get(sp, "kspace_other")
        print(f"  {sp:<22s} [{cls:<14s}] {v['wall_sec']:10.4f} s  "
              f"{v['pct_of_kspace_wall']:6.2f}% of kspace "
              f"({v['events']} events, {v['bytes']/1e6:.2f} MB)")
    print("\n--- Summary ---")
    print(f"  Kspace comm / total comm wall ......... "
          f"{s['kspace_pct_of_comm_wall']:6.2f}%")
    print(f"  Kspace GLOBAL comm / kspace comm wall . "
          f"{s['global_pct_of_kspace_wall']:6.2f}%   "
          f"({s['kspace_global_wall_sec']:.4f} s, "
          f"{s['kspace_global_bytes']/1e6:.2f} MB)")
    print(f"  Kspace LOCAL  comm / kspace comm wall . "
          f"{s['local_pct_of_kspace_wall']:6.2f}%   "
          f"({s['kspace_local_wall_sec']:.4f} s, "
          f"{s['kspace_local_bytes']/1e6:.2f} MB)")
    print(f"  GLOBAL comm / total comm wall ......... "
          f"{s['global_pct_of_comm_wall']:6.2f}%")
    print()


def trend_row(r):
    s = r["summary"]
    return {
        "ranks": r["num_ranks"],
        "comm_wall_sec": round(r["total_comm"]["wall_sec"], 6),
        "kspace_pct_of_comm": round(s["kspace_pct_of_comm_wall"], 2),
        "kspace_global_pct_of_kspace": round(s["global_pct_of_kspace_wall"], 2),
        "kspace_local_pct_of_kspace": round(s["local_pct_of_kspace_wall"], 2),
        "global_pct_of_comm": round(s["global_pct_of_comm_wall"], 2),
        "kspace_global_wall_sec": round(s["kspace_global_wall_sec"], 6),
        "kspace_local_wall_sec": round(s["kspace_local_wall_sec"], 6),
        "kspace_global_events": s["kspace_global_events"],
        "kspace_local_events": s["kspace_local_events"],
    }


def find_run_files(run_dir):
    ccdgs = sorted(glob.glob(os.path.join(run_dir, "*.ccdg")))
    phase = os.path.join(run_dir, "phase_trace.csv")
    if not ccdgs:
        return None
    if not os.path.isfile(phase):
        return None
    return ccdgs[0], phase


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Marker-driven global/local comm attribution "
                    "(phase trace x CCDG wall time)")
    ap.add_argument("ccdg", nargs="?", help="CCDG JSON file (single mode)")
    ap.add_argument("phase_csv", nargs="?",
                    help="phase_trace.csv (single mode)")
    ap.add_argument("--batch", nargs="+", metavar="RUN_DIR",
                    help="run directories, each containing *.ccdg and "
                         "phase_trace.csv")
    ap.add_argument("--json", help="write analysis JSON to this path")
    ap.add_argument("--csv", help="write strong-scaling trend CSV "
                                   "(batch mode)")
    args = ap.parse_args()

    if args.batch:
        rows = []
        results = []
        for d in args.batch:
            found = find_run_files(d)
            if not found:
                print(f"WARNING: skip {d} (missing *.ccdg or phase_trace.csv)",
                      file=sys.stderr)
                continue
            r = analyze(found[0], found[1])
            results.append(r)
            rows.append(trend_row(r))
            print_report(r)
        rows.sort(key=lambda x: x["ranks"])
        if rows:
            print("=" * 72)
            print("Strong-scaling trend (fixed problem size)")
            print("=" * 72)
            hdr = ("ranks", "kspace%", "kspace-global%of-kspace",
                   "global%of-comm", "global_wall_s", "local_wall_s")
            print(f"  {hdr[0]:>6s} {hdr[1]:>9s} {hdr[2]:>22s} "
                  f"{hdr[3]:>15s} {hdr[4]:>13s} {hdr[5]:>13s}")
            for x in rows:
                print(f"  {x['ranks']:>6d} "
                      f"{x['kspace_pct_of_comm']:>8.2f}% "
                      f"{x['kspace_global_pct_of_kspace']:>21.2f}% "
                      f"{x['global_pct_of_comm']:>14.2f}% "
                      f"{x['kspace_global_wall_sec']:>13.4f} "
                      f"{x['kspace_local_wall_sec']:>13.4f}")
            if args.csv:
                with open(args.csv, "w", newline="") as f:
                    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                    w.writeheader()
                    w.writerows(rows)
                print(f"\nTrend CSV written to {args.csv}")
        if args.json:
            with open(args.json, "w") as f:
                json.dump(results, f, indent=2)
            print(f"Analysis JSON written to {args.json}")
    elif args.ccdg and args.phase_csv:
        r = analyze(args.ccdg, args.phase_csv)
        print_report(r)
        if args.json:
            with open(args.json, "w") as f:
                json.dump(r, f, indent=2)
            print(f"Analysis JSON written to {args.json}")
    else:
        ap.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
