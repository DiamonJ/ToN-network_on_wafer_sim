#!/usr/bin/env python3
"""Large-scale BookSim synthetic traffic and WSE gating experiments.

Produces three reproducible suites:

* uniform_scale: 4x4 through 256x256, injection rate vs queue/E2E latency
* gating_scale: WSE phase width vs queue/E2E latency across the same meshes
* traffic_compare: 4x4 uniform, transpose, and uniform+WSE vs injection rate

Large meshes use a bounded sampling protocol: one warm-up period, one measured
period, and a final drain.  Every point records its sampling period and status;
SATURATED/TIMEOUT points are retained instead of truncating the whole curve.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import tempfile
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional, Tuple


ROOT = Path(__file__).resolve().parent
BOOKSIM = ROOT / "src" / "booksim"
CFG = ROOT / "traffic_wse_4x4.cfg"
DEFAULT_OUTPUT = ROOT / "sweep_results" / "synthetic_scaling"
MESH_SIDES = [4, 8, 16, 32, 64, 128, 256]
# packet injection_rate * mesh_side; keeps offered bisection load comparable.
NORMALIZED_LOADS = [0.0256, 0.0512, 0.1024, 0.1536, 0.2048, 0.3072, 0.3840]
COMMON_RATES = [0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.0012, 0.0015]
PHASE_WIDTHS = [2, 8, 32, 64, 128, 512]
COMPARE_RATES = [round(0.01 + 0.005 * i, 3) for i in range(35)]  # 0.010..0.180


def finite(value: float) -> bool:
    return math.isfinite(value)


def sample_period(k: int, injection_rate: float, target_packets: int = 1700) -> int:
    """Bound work while retaining enough injected traffic for useful stats."""
    expected_per_cycle = max(1e-12, k * k * injection_rate)
    needed = math.ceil(target_packets / expected_per_cycle)
    minimum = 64 if k >= 256 else 512
    scale_cap = k if k >= 256 else (4 * k if k >= 128 else 20000)
    return max(minimum, min(20000, needed, scale_cap))


def point_id(spec: dict[str, Any]) -> str:
    fields = [
        spec["suite"], f"k{spec['k']}", spec["traffic"],
        f"inj{spec['injection_rate']:.10g}",
        f"g{spec['wse_gating']}", f"pw{spec['phase_width']}",
    ]
    return "_".join(fields).replace(".", "p")


def parse_result(stdout: str) -> Optional[Tuple[float, float, float, float]]:
    lines = [line for line in stdout.splitlines() if line.startswith("results:")]
    if not lines:
        return None
    parts = lines[-1].split(",")
    if len(parts) < 9:
        return None
    try:
        avg_plat = float(parts[5])
        avg_nlat = float(parts[8])
        reported_packet_injection = float(parts[3])
    except ValueError:
        return None
    if not finite(avg_plat) or not finite(avg_nlat):
        return None
    return max(0.0, avg_plat - avg_nlat), avg_plat, avg_nlat, reported_packet_injection


def run_point(spec: dict[str, Any], output: Path, timeout: int,
              force: bool = False, retry_failures: bool = False) -> dict[str, Any]:
    pid = point_id(spec)
    cache = output / "points" / f"{pid}.json"
    log_path = output / "logs" / f"{pid}.log"
    if cache.exists() and not force:
        cached = json.loads(cache.read_text(encoding="utf-8"))
        if not retry_failures or cached.get("status") == "PASS":
            return cached

    legacy_compare = spec["suite"] == "traffic_compare"
    period = 1000 if legacy_compare else sample_period(
        spec["k"], spec["injection_rate"],
        target_packets=(512 if spec["suite"] == "uniform_common" else 1700),
    )
    if spec["suite"] == "gating_scale" and spec["k"] >= 256:
        period = 64
    # Four VCs keep 65,536-router construction tractable and still provide
    # one VC per WSE direction class. The 4x4 legacy comparison overrides it.
    num_vcs = 24 if legacy_compare else 4
    vc_buf = 8 if legacy_compare else 4
    warmup_periods = 3 if legacy_compare else 1
    max_samples = 50 if legacy_compare else 2
    latency_threshold = 10000.0 if legacy_compare else 1000000000.0
    command = [
        str(BOOKSIM), str(CFG),
        f"k={spec['k']}", "n=2", f"traffic={spec['traffic']}",
        f"injection_rate={spec['injection_rate']:.12g}",
        # Preserve the original sweep_queue_time.py semantics: injection_rate
        # is packets/node/cycle (BookSim's default), not flits/node/cycle.
        "injection_rate_uses_flits=0", "packet_size=8",
        f"num_vcs={num_vcs}", f"vc_buf_size={vc_buf}",
        f"wse_gating={spec['wse_gating']}",
        f"wse_phase_width={spec['phase_width']}",
        "wse_strip_width=2", "wse_vc_split=1",
        "include_queuing=1", "sim_count=1", "seed=0",
        f"warmup_periods={warmup_periods}", f"sample_period={period}",
        f"max_samples={max_samples}", f"latency_thres={latency_threshold}",
        "print_csv_results=1",
        "stats_out=/dev/null",
    ]
    started = time.monotonic()
    timed_out = False
    try:
        proc = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=timeout)
        stdout, rc = proc.stdout, proc.returncode
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        rc = 124
        if isinstance(exc.stdout, bytes):
            stdout = exc.stdout.decode(errors="replace")
        else:
            stdout = exc.stdout or ""
    elapsed = time.monotonic() - started
    log_path.write_text(stdout, encoding="utf-8")
    parsed = parse_result(stdout)
    if timed_out:
        status = "TIMEOUT"
    elif parsed is None and ("Simulation unstable" in stdout or
                             "exceeded" in stdout or
                             "Too many sample periods" in stdout):
        status = "SATURATED"
    elif parsed is None:
        status = "NO_SAMPLES"
    else:
        status = "PASS"
    queue, e2e, network, reported = parsed if parsed else (None, None, None, None)
    row = {
        **spec,
        "nodes": spec["k"] * spec["k"],
        "mesh": f"{spec['k']}x{spec['k']}",
        "normalized_load": spec["injection_rate"] * spec["k"],
        "full_wavelet_period_cycles": (
            4 * spec["phase_width"] if spec["wse_gating"] else 0
        ),
        "avg_queue_cycles": queue,
        "avg_end_to_end_cycles": e2e,
        "avg_network_cycles": network,
        "reported_packet_injection_rate": reported,
        "sample_period": period,
        "warmup_periods": warmup_periods,
        "max_samples": max_samples,
        "sampling_protocol": "converged" if legacy_compare else "bounded_one_period",
        "num_vcs": num_vcs,
        "vc_buf_size": vc_buf,
        "packet_size_flits": 8,
        "status": status,
        "return_code": rc,
        "host_elapsed_seconds": elapsed,
        "log": str(log_path),
    }
    cache.write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    print(
        f"{spec['suite']:<15} k={spec['k']:3d} inj={spec['injection_rate']:.6g} "
        f"traffic={spec['traffic']:<9} pw={spec['phase_width']:3d} "
        f"{status:<10} queue={queue} e2e={e2e} host={elapsed:.1f}s",
        flush=True,
    )
    return row


def specs_for(suites: Iterable[str], max_k: int) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    meshes = [k for k in MESH_SIDES if k <= max_k]
    if "uniform_scale" in suites:
        for k in meshes:
            for load in NORMALIZED_LOADS:
                specs.append({
                    "suite": "uniform_scale", "k": k, "traffic": "uniform",
                    "injection_rate": load / k, "wse_gating": 0, "phase_width": 64,
                })
    if "uniform_common" in suites:
        for k in meshes:
            for rate in COMMON_RATES:
                specs.append({
                    "suite": "uniform_common", "k": k, "traffic": "uniform",
                    "injection_rate": rate, "wse_gating": 0, "phase_width": 64,
                })
    if "gating_scale" in suites:
        for k in meshes:
            for width in PHASE_WIDTHS:
                specs.append({
                    "suite": "gating_scale", "k": k, "traffic": "uniform",
                    # 256x256 reaches the practical drain boundary at 0.0004;
                    # use the highest empirically completed rate for its
                    # gating sensitivity curve.
                    "injection_rate": (0.0001 if k >= 256 else 0.1024 / k),
                    "wse_gating": 1,
                    "phase_width": width,
                })
    if "traffic_compare" in suites:
        for rate in COMPARE_RATES:
            scenarios = []
            if rate <= 0.180:
                scenarios.append((4, "uniform", 0))
            if rate <= 0.090:
                scenarios.append((8, "uniform", 0))
            if rate <= 0.085:
                scenarios.append((4, "transpose", 0))
            if rate <= 0.130:
                scenarios.append((4, "uniform", 1))
            for k, traffic, gating in scenarios:
                specs.append({
                    "suite": "traffic_compare", "k": k, "traffic": traffic,
                    "injection_rate": rate, "wse_gating": gating, "phase_width": 64,
                })
    return specs


def plot_lines(rows: list[dict[str, Any]], output: Path, suite: str,
               xkey: str, ykey: str, xlabel: str, ylabel: str,
               filename: str, group_label, max_x: Optional[float] = None,
               explicit_xticks: Optional[list] = None,
               y_scale: Optional[str] = None,
               saturation_by_group: Optional[dict] = None) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter, NullFormatter
    except ImportError:
        return
    selected = [r for r in rows if r["suite"] == suite and
                r["status"] in ("PASS", "ESTIMATED")]
    if max_x is not None:
        selected = [r for r in selected if float(r[xkey]) <= max_x]
    if not selected:
        return
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in selected:
        groups[group_label(row)].append(row)
    fig, ax = plt.subplots(figsize=(10, 6.2))
    markers = ["o", "s", "^", "D", "v", "P", "h"]
    def group_order(item):
        label, group = item
        return (int(group[0].get("k", 0)), label)

    for idx, (label, group) in enumerate(sorted(groups.items(), key=group_order)):
        group.sort(key=lambda row: row[xkey])
        color = f"C{idx % 10}"
        saturation = (saturation_by_group or {}).get(label)
        if saturation is not None:
            group = [row for row in group if float(row[xkey]) < saturation]
            ax.axvline(saturation, color=color, linestyle="--", linewidth=1.0,
                       alpha=0.45)
            display_label = f"{label} (sat≈{saturation:g})"
        else:
            display_label = label
        line_style = "-"
        if suite == "traffic_compare":
            if label.endswith("transpose"):
                line_style = "--"
            elif label.startswith("uniform + WSE"):
                line_style = "-."
        ax.plot([r[xkey] for r in group], [r[ykey] for r in group],
                linewidth=1.8, linestyle=line_style, color=color,
                label=display_label)
        measured = [r for r in group if r["status"] == "PASS"]
        estimated = [r for r in group if r["status"] == "ESTIMATED"]
        if suite == "traffic_compare" and label.endswith("transpose"):
            # Uniform and transpose source-queue values are intentionally
            # identical below congestion. A large hollow marker exposes the
            # filled uniform square underneath without changing either value.
            ax.scatter([r[xkey] for r in measured], [r[ykey] for r in measured],
                       marker="o", s=70, facecolors="none", edgecolors=color,
                       linewidths=1.8, zorder=5)
        else:
            ax.scatter([r[xkey] for r in measured], [r[ykey] for r in measured],
                       marker=markers[idx % len(markers)], s=32,
                       color=color, zorder=3)
        if estimated:
            ax.scatter([r[xkey] for r in estimated], [r[ykey] for r in estimated],
                       marker="X", s=55, color=color, edgecolors="black",
                       linewidths=0.5, zorder=4)
    if xkey == "injection_rate":
        ax.set_xscale("log")
        if explicit_xticks:
            ax.set_xticks(explicit_xticks)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
    elif xkey == "phase_width":
        ax.set_xscale("log", base=2)
        ax.set_xticks(PHASE_WIDTHS)
        ax.set_xticklabels([str(value) for value in PHASE_WIDTHS])
        ax.xaxis.set_minor_formatter(NullFormatter())
    elif explicit_xticks:
        ax.set_xticks(explicit_xticks)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    if y_scale == "linear":
        ax.set_ylim(bottom=0)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    elif y_scale == "log":
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
    elif ykey == "avg_queue_cycles":
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9, ncol=2)
    if any(r["status"] == "ESTIMATED" for r in selected):
        fig.text(0.5, 0.005, "X markers: extrapolated from 128x128 and calibrated by 256x256 BookSim anchors",
                 ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(output / f"{filename}.png", dpi=170)
    fig.savefig(output / f"{filename}.svg")
    plt.close(fig)


def make_plots(rows: list[dict[str, Any]], output: Path) -> None:
    mesh_label = lambda r: r["mesh"]
    # Preserve the absolute-rate, common stable-window figures as appendices.
    plot_lines(rows, output, "uniform_common", "injection_rate", "avg_queue_cycles",
               "Injection rate (packets/node/cycle)", "Average queueing delay (cycles)",
               "uniform_mesh_common_injection_vs_queue", mesh_label,
               explicit_xticks=COMMON_RATES, y_scale="linear")
    plot_lines(rows, output, "uniform_common", "injection_rate", "avg_end_to_end_cycles",
               "Injection rate (packets/node/cycle)", "Average packet end-to-end latency (cycles)",
               "uniform_mesh_common_injection_vs_e2e", mesh_label,
               explicit_xticks=COMMON_RATES, y_scale="log")
    saturation_by_mesh = {
        f"{k}x{k}": min(0.125, 0.4 / k) for k in MESH_SIDES
    }
    # Primary strict steady-state figures: raw injection rate, with each curve
    # ending before its topology-specific saturation rate. Dashed vertical
    # lines mark saturation; no finite latency point is drawn at/after it.
    plot_lines(rows, output, "uniform_scale", "injection_rate", "avg_queue_cycles",
               "Injection rate (packets/node/cycle)",
               "Average queueing delay (cycles)",
               "uniform_mesh_injection_vs_queue", mesh_label,
               y_scale=None, saturation_by_group=saturation_by_mesh)
    plot_lines(rows, output, "uniform_scale", "injection_rate", "avg_end_to_end_cycles",
               "Injection rate (packets/node/cycle)",
               "Average packet end-to-end latency (cycles)",
               "uniform_mesh_injection_vs_e2e", mesh_label,
               y_scale="log", saturation_by_group=saturation_by_mesh)
    # Capacity-normalized versions are retained as appendices.
    plot_lines(rows, output, "uniform_scale", "normalized_load", "avg_queue_cycles",
               "Normalized injection load (packet rate × mesh side)",
               "Average queueing delay (cycles)",
               "uniform_mesh_normalized_injection_vs_queue", mesh_label,
               explicit_xticks=NORMALIZED_LOADS)
    plot_lines(rows, output, "uniform_scale", "normalized_load", "avg_end_to_end_cycles",
               "Normalized injection load (packet rate × mesh side)",
               "Average packet end-to-end latency (cycles)",
               "uniform_mesh_normalized_injection_vs_e2e", mesh_label,
               explicit_xticks=NORMALIZED_LOADS, y_scale="log")
    plot_lines(rows, output, "gating_scale", "phase_width", "avg_queue_cycles",
               "WSE phase width (cycles)", "Average queueing delay (cycles)",
               "wse_gating_vs_queue", mesh_label)
    plot_lines(rows, output, "gating_scale", "phase_width", "avg_end_to_end_cycles",
               "WSE phase width (cycles)", "Average packet end-to-end latency (cycles)",
               "wse_gating_vs_e2e", mesh_label)

    def compare_label(row: dict[str, Any]) -> str:
        if row["wse_gating"]:
            return f"{row['mesh']} uniform + WSE gating (pw=64)"
        return f"{row['mesh']} {row['traffic']}"
    compare_saturation = {
        "4x4 uniform": 0.185,
        "8x8 uniform": 0.095,
        "4x4 transpose": 0.090,
        "4x4 uniform + WSE gating (pw=64)": 0.135,
    }
    compare_ticks = [0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.18]
    plot_lines(rows, output, "traffic_compare", "injection_rate", "avg_queue_cycles",
               "Injection rate (packets/node/cycle)", "Average queueing delay (cycles)",
               "traffic_compare_injection_vs_queue", compare_label,
               explicit_xticks=compare_ticks,
               saturation_by_group=compare_saturation)
    plot_lines(rows, output, "traffic_compare", "injection_rate", "avg_end_to_end_cycles",
               "Injection rate (packets/node/cycle)", "Average packet end-to-end latency (cycles)",
               "traffic_compare_injection_vs_e2e", compare_label,
               explicit_xticks=compare_ticks,
               saturation_by_group=compare_saturation)


def add_large_mesh_estimates(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Fill expensive 256x256 points from 128x128 plus measured anchors.

    The estimator preserves measured source-queue delay at the same normalized
    offered load. Network-latency increments are scaled by the path-length
    ratio 256/128, then anchored to the measured 256x256 low-load latency.
    """
    for row in rows:
        row.setdefault("provenance", "booksim_cycle_accurate")
        row.setdefault("execution_status", row["status"])

    def key_load(row: dict[str, Any]) -> float:
        return round(float(row["normalized_load"]), 7)

    uniform_128 = {
        key_load(row): row for row in rows
        if row["suite"] == "uniform_scale" and row["k"] == 128
        and row["status"] == "PASS"
    }
    uniform_256 = [row for row in rows
                   if row["suite"] == "uniform_scale" and row["k"] == 256]
    anchors_128 = sorted(uniform_128.values(), key=lambda row: row["normalized_load"])
    anchors_256 = sorted((row for row in uniform_256 if row["status"] == "PASS"),
                         key=lambda row: row["normalized_load"])
    if anchors_128 and anchors_256:
        source_network0 = float(anchors_128[0]["avg_network_cycles"])
        target_network0 = float(anchors_256[0]["avg_network_cycles"])
        for row in uniform_256:
            if row["status"] == "PASS":
                continue
            source = uniform_128.get(key_load(row))
            if not source:
                continue
            row["execution_status"] = row["status"]
            row["status"] = "ESTIMATED"
            row["provenance"] = "128x128_scaled_with_256x256_anchor"
            row["avg_queue_cycles"] = source["avg_queue_cycles"]
            row["avg_network_cycles"] = (
                target_network0 + 2.0 *
                (float(source["avg_network_cycles"]) - source_network0)
            )
            row["avg_end_to_end_cycles"] = (
                row["avg_network_cycles"] + row["avg_queue_cycles"]
            )
        # Finite-window anchors carry about 1% sampling noise. Preserve the
        # physically monotone injection-latency relationship for extrapolated
        # points without changing any measured PASS value.
        running_e2e = 0.0
        for row in sorted(uniform_256, key=lambda item: item["injection_rate"]):
            value = float(row["avg_end_to_end_cycles"])
            if row["status"] == "ESTIMATED" and value < running_e2e:
                row["avg_end_to_end_cycles"] = running_e2e
                row["avg_network_cycles"] = running_e2e - float(row["avg_queue_cycles"])
            running_e2e = max(running_e2e, float(row["avg_end_to_end_cycles"]))

    gating_128 = {
        int(row["phase_width"]): row for row in rows
        if row["suite"] == "gating_scale" and row["k"] == 128
        and row["status"] == "PASS"
    }
    gating_256 = {
        int(row["phase_width"]): row for row in rows
        if row["suite"] == "gating_scale" and row["k"] == 256
    }
    target_network0 = (float(anchors_256[0]["avg_network_cycles"])
                       if anchors_256 else 854.5)
    for width in PHASE_WIDTHS:
        if width in gating_256 and gating_256[width]["status"] == "PASS":
            continue
        source = gating_128.get(width)
        if not source:
            continue
        if width in gating_256:
            row = gating_256[width]
            row["execution_status"] = row["status"]
        else:
            row = dict(source)
            rows.append(row)
            row.update({
                "k": 256, "nodes": 65536, "mesh": "256x256",
                "injection_rate": 0.0001, "normalized_load": 0.0256,
                "sample_period": 64, "host_elapsed_seconds": None,
                "log": "", "return_code": None,
                "execution_status": "NOT_RUN",
            })
        row["status"] = "ESTIMATED"
        row["provenance"] = "128x128_gating_scaled_with_256x256_anchor"
        row["avg_queue_cycles"] = source["avg_queue_cycles"]
        row["avg_network_cycles"] = target_network0
        row["avg_end_to_end_cycles"] = target_network0 + float(source["avg_queue_cycles"])
    return rows


def add_common_rate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Complete the common-x uniform curves for 128x128 and 256x256."""
    existing = {(int(r["k"]), round(float(r["injection_rate"]), 7)): r
                for r in rows if r["suite"] == "uniform_common"}
    source_by_k: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["suite"] == "uniform_scale" and row["status"] in ("PASS", "ESTIMATED"):
            source_by_k[int(row["k"])].append(row)
    for source in source_by_k.values():
        source.sort(key=lambda row: float(row["injection_rate"]))

    for k in MESH_SIDES:
        source = source_by_k.get(k, [])
        if not source:
            continue
        for rate in COMMON_RATES:
            key = (k, round(rate, 7))
            if key in existing:
                continue
            exact = next((row for row in source
                          if math.isclose(float(row["injection_rate"]), rate,
                                          rel_tol=0.0, abs_tol=1e-12)), None)
            if exact:
                row = dict(exact)
                row["suite"] = "uniform_common"
                rows.append(row)
                existing[key] = row
                continue
            left = max((row for row in source if float(row["injection_rate"]) < rate),
                       key=lambda row: float(row["injection_rate"]), default=None)
            right = min((row for row in source if float(row["injection_rate"]) > rate),
                        key=lambda row: float(row["injection_rate"]), default=None)
            template = left or right
            if template is None:
                continue
            row = dict(template)
            if left and right:
                x0, x1 = float(left["injection_rate"]), float(right["injection_rate"])
                weight = (rate - x0) / (x1 - x0)
                def blend(field: str) -> float:
                    return float(left[field]) + weight * (float(right[field]) - float(left[field]))
                queue = blend("avg_queue_cycles")
                network = blend("avg_network_cycles")
            elif right:
                # Below the first measured load: use zero-load queueing and
                # the first measured propagation latency.
                queue = 0.0
                network = float(right["avg_network_cycles"])
            else:
                queue = float(left["avg_queue_cycles"])
                network = float(left["avg_network_cycles"])
            row.update({
                "suite": "uniform_common", "injection_rate": rate,
                "normalized_load": rate * k,
                "avg_queue_cycles": queue,
                "avg_network_cycles": network,
                "avg_end_to_end_cycles": network + queue,
                "status": "ESTIMATED",
                "execution_status": "INTERPOLATED",
                "provenance": "interpolated_from_uniform_scale",
                "host_elapsed_seconds": None, "log": "", "return_code": None,
            })
            rows.append(row)
            existing[key] = row
    return rows


def write_outputs(rows: list[dict[str, Any]], output: Path) -> None:
    rows = add_large_mesh_estimates(rows)
    rows = add_common_rate_rows(rows)
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with (output / "all_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    (output / "all_results.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    for suite in ("uniform_common", "uniform_scale", "gating_scale", "traffic_compare"):
        subset = [r for r in rows if r["suite"] == suite]
        if not subset:
            continue
        with (output / f"{suite}.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(subset)
    status_counts = defaultdict(int)
    for row in rows:
        status_counts[row["status"]] += 1
    readme = [
        "# Synthetic scaling experiment results", "",
        "Metrics: queue = average packet latency - average network latency; end-to-end = average packet latency.",
        "Injection rate is in packets/node/cycle, matching the original sweep_queue_time.py semantics. Large-scale runs use one warm-up and one measured period followed by complete drain.",
        "4x4 traffic comparison uses 24 VCs x 8 flits; mesh/gating scaling uses 4 VCs x 4 flits for 65,536-router feasibility.",
        "", "Status counts: " + ", ".join(f"{k}={v}" for k, v in sorted(status_counts.items())),
        "", "256x256 X markers are extrapolated, not cycle-accurate points: uniform uses matching 128x128 normalized-load queueing plus a path-length-scaled network increment anchored by measured 256x256 latency; missing gating points use the 128x128 gating penalty plus the measured 256x256 network baseline.",
        "", "See all_results.csv for per-point sample_period, host runtime, and log path.", "",
    ]
    (output / "README.md").write_text("\n".join(readme), encoding="utf-8")
    make_plots(rows, output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=["all", "uniform_common", "uniform_scale", "gating_scale", "traffic_compare"], default="all")
    parser.add_argument("--max-k", type=int, default=256)
    parser.add_argument("--min-k", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--retry-failures", action="store_true")
    parser.add_argument("--plot-only", action="store_true",
                        help="skip simulation and rebuild tables/plots from point cache")
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    # Keep Matplotlib cache portable and writable on Linux/macOS/Windows.
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "dpmd-matplotlib-cache")
    )
    suites = (["uniform_common", "uniform_scale", "gating_scale", "traffic_compare"]
              if args.suite == "all" else [args.suite])
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "points").mkdir(exist_ok=True)
    (args.output / "logs").mkdir(exist_ok=True)
    if not args.plot_only:
        specs = [spec for spec in specs_for(suites, args.max_k)
                 if spec["k"] >= args.min_k]
        rows = [run_point(spec, args.output, args.timeout, args.force,
                          args.retry_failures) for spec in specs]
    # Merge cached suites so separate invocations regenerate the complete plots.
    cached = {}
    for path in (args.output / "points").glob("*.json"):
        row = json.loads(path.read_text(encoding="utf-8"))
        cached[point_id(row)] = row
    if args.plot_only and not cached:
        saved = args.output / "all_results.json"
        if not saved.is_file():
            raise SystemExit(f"plot cache/data not found: {saved}")
        for row in json.loads(saved.read_text(encoding="utf-8")):
            cached[point_id(row)] = row
    rows = sorted(cached.values(), key=lambda r: (
        r["suite"], r["k"], r["traffic"], r["wse_gating"],
        r["injection_rate"], r["phase_width"],
    ))
    write_outputs(rows, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
