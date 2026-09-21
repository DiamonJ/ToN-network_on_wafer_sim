#!/usr/bin/env python3
"""Run and summarize the Phase-3 WSE experiment matrix."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
DEFAULT_CASES = {
    16: ROOT / "runs/pipeline/short_lialocl_2688a_16r_20260910_093912",
    64: ROOT / "cases/wse_short_64",
    256: ROOT / "cases/wse_short_256",
}
BASELINE_16 = DEFAULT_CASES[16] / "evaluation.json"
HISTORICAL_CSV = ROOT / "booksim2/results/ccdg_mesh_results.csv"


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def run(command: list[str], log_path: Path) -> None:
    result = subprocess.run(
        command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    log_path.write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}); see {log_path}")


def parse_stat(text: str, name: str, cast: type = int) -> Any:
    match = re.search(rf"^{re.escape(name)} = ([0-9.eE+-]+);", text, re.MULTILINE)
    if not match:
        raise ValueError(f"missing BookSim stat {name}")
    return cast(match.group(1))


def flatten_stages(
    program: dict[str, Any], replay: dict[str, Any], replay_est: Path, stats_text: str
) -> list[dict[str, Any]]:
    release_by_node = {}
    for line in replay_est.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            node, release = line.split()
            release_by_node[int(node)] = int(release)
    stage_releases: dict[int, list[int]] = defaultdict(list)
    for node in replay["nodes"]:
        if node.get("wse_kind") == "branch":
            stage_releases[int(node["wse_stage_idx"])].append(
                release_by_node[int(node["id"])]
            )

    rows = []
    stage_idx = 0
    for timestep in program["timesteps"]:
        for stage in timestep["stages"]:
            waves = stage["wavefronts"]
            branches = sum(
                1
                for wave in waves
                for slot in wave["footprint_slots"]
                if re.fullmatch(r"\d+,\d+->\d+,\d+", slot["link"])
            )
            link_bytes = sum(
                int(slot["bytes"])
                for wave in waves
                for slot in wave["footprint_slots"]
                if re.fullmatch(r"\d+,\d+->\d+,\d+", slot["link"])
            )
            expected = parse_stat(
                stats_text, f"wse_stage_{stage_idx}_branches_expected"
            )
            delivered = parse_stat(
                stats_text, f"wse_stage_{stage_idx}_branches_delivered"
            )
            completion = parse_stat(
                stats_text, f"wse_stage_{stage_idx}_completion_cycle"
            )
            release = min(stage_releases[stage_idx])
            rows.append(
                {
                    "stage_index": stage_idx,
                    "stage_id": stage["id"],
                    "phase": stage["phase"],
                    "axis": stage["axis"],
                    "mode": stage["mode"],
                    "b": stage["b"],
                    "wavefronts": len(waves),
                    "branches": branches,
                    "branches_expected": expected,
                    "branches_delivered": delivered,
                    "link_bytes": link_bytes,
                    "compiler_start_cycle": stage["start_cycle"],
                    "compiler_completion_cycle": stage["completion_cycle"],
                    "booksim_release_cycle": release,
                    "booksim_completion_cycle": completion,
                    "booksim_stage_span_cycles": completion - release + 1,
                }
            )
            stage_idx += 1
    return rows


def link_rows(
    program: dict[str, Any], rank_count: int, profile: str, fold: str
) -> list[dict[str, Any]]:
    loads: dict[tuple[int, int, int, int], int] = defaultdict(int)
    for timestep in program["timesteps"]:
        for stage in timestep["stages"]:
            for wave in stage["wavefronts"]:
                for slot in wave["footprint_slots"]:
                    match = re.fullmatch(r"(\d+),(\d+)->(\d+),(\d+)", slot["link"])
                    if match:
                        loads[tuple(map(int, match.groups()))] += int(slot["flits"])
    rows = []
    for (sx, sy, dx, dy), flits in sorted(loads.items()):
        rows.append(
            {
                "ranks": rank_count,
                "profile": profile,
                "fold": fold,
                "src_x": sx,
                "src_y": sy,
                "dst_x": dx,
                "dst_y": dy,
                "direction": (
                    "H+" if dx > sx else "H-" if dx < sx else "V+" if dy > sy else "V-"
                ),
                "flit_cycles": flits,
            }
        )
    return rows


def historical_baselines() -> list[dict[str, Any]]:
    rows = []
    if BASELINE_16.exists():
        data = read_json(BASELINE_16)
        for name, result in data["sims"].items():
            rows.append(
                {
                    "ranks": 16,
                    "name": name,
                    "cycles": int(result["cycles_per_iter"]),
                    "steps_per_sec": float(result["timesteps_per_sec"]),
                    "source": str(BASELINE_16),
                }
            )
    if HISTORICAL_CSV.exists():
        with HISTORICAL_CSV.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if (
                    row.get("mode") == "short"
                    and row.get("ranks") in ("64", "256")
                    and "lialocl_2688a" in row.get("ccdg_dir", "")
                    and "cap2.5e10" in row.get("ccdg_dir", "")
                    and row.get("gating") == "free"
                ):
                    rows.append(
                        {
                            "ranks": int(row["ranks"]),
                            "name": "free_historical",
                            "cycles": int(float(row["cycles_per_iter"])),
                            "steps_per_sec": float(row["timesteps_per_sec"]),
                            "source": str(HISTORICAL_CSV),
                        }
                    )
    # Keep the newest matching historical row for each scale/name.
    dedup = {(row["ranks"], row["name"]): row for row in rows}
    return sorted(dedup.values(), key=lambda row: (row["ranks"], row["name"]))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def render_heatmaps(output: Path, results: list[dict[str, Any]], links: list[dict[str, Any]]) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return
    grouped: dict[tuple[int, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in links:
        grouped[(row["ranks"], row["profile"], row["fold"])].append(row)
    result_by_key = {
        (row["ranks"], row["profile"], row["fold"]): row for row in results
    }
    for key, group in grouped.items():
        ranks, profile, fold = key
        k = round(ranks**0.5)
        horizontal = np.zeros((k, k))
        vertical = np.zeros((k, k))
        window = max(
            1,
            result_by_key[key]["compiler_cycles"]
            - result_by_key[key]["compute_barrier_cycles"],
        )
        for row in group:
            matrix = horizontal if row["direction"].startswith("H") else vertical
            matrix[row["src_y"], row["src_x"]] += row["flit_cycles"] / window
        vmax = max(float(horizontal.max()), float(vertical.max()), 1e-12)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        for axis, matrix, title in zip(
            axes, (horizontal, vertical), ("Horizontal directed links", "Vertical directed links")
        ):
            image = axis.imshow(matrix, origin="lower", vmin=0, vmax=vmax, cmap="viridis")
            axis.set_title(title)
            axis.set_xlabel("Mesh X coordinate")
            axis.set_ylabel("Mesh Y coordinate")
        fig.colorbar(image, ax=axes, label="Aggregate flit cycles / communication window")
        fig.suptitle(f"{ranks} ranks · {profile} · fold {fold}")
        fig.savefig(output / f"heatmap_r{ranks}_{profile}_fold_{fold}.png", dpi=160)
        plt.close(fig)


def render_svg_heatmaps(
    output: Path, results: list[dict[str, Any]], links: list[dict[str, Any]]
) -> None:
    grouped: dict[tuple[int, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in links:
        grouped[(row["ranks"], row["profile"], row["fold"])].append(row)
    result_by_key = {
        (row["ranks"], row["profile"], row["fold"]): row for row in results
    }
    for key, group in grouped.items():
        ranks, profile, fold = key
        k = round(ranks**0.5)
        window = max(
            1,
            result_by_key[key]["compiler_cycles"]
            - result_by_key[key]["compute_barrier_cycles"],
        )
        matrices = {"Horizontal links": defaultdict(float), "Vertical links": defaultdict(float)}
        for row in group:
            panel = "Horizontal links" if row["direction"].startswith("H") else "Vertical links"
            matrices[panel][(row["src_x"], row["src_y"])] += row["flit_cycles"] / window
        vmax = max((value for matrix in matrices.values() for value in matrix.values()), default=1.0)
        cell = max(8, min(24, 320 // k))
        panel_width = k * cell
        width, height = panel_width * 2 + 120, k * cell + 100
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">',
            '<rect width="100%" height="100%" fill="white"/>',
            f'<text x="20" y="24" font-family="sans-serif" font-size="16">'
            f'{ranks} ranks · {profile} · fold {fold}</text>',
            '<text x="20" y="44" font-family="sans-serif" font-size="11">'
            'Cell color: aggregate directed-link flit cycles / communication window</text>',
        ]
        for panel_idx, (title, matrix) in enumerate(matrices.items()):
            origin_x = 40 + panel_idx * (panel_width + 60)
            origin_y = 72
            parts.append(
                f'<text x="{origin_x}" y="{origin_y - 8}" font-family="sans-serif" '
                f'font-size="12">{title}</text>'
            )
            for y in range(k):
                for x in range(k):
                    value = matrix[(x, y)]
                    opacity = 0.05 + 0.95 * (value / vmax if vmax else 0.0)
                    draw_y = origin_y + (k - 1 - y) * cell
                    parts.append(
                        f'<rect x="{origin_x + x * cell}" y="{draw_y}" '
                        f'width="{cell}" height="{cell}" fill="rgb(33,102,172)" '
                        f'fill-opacity="{opacity:.4f}" stroke="rgb(220,220,220)" '
                        f'stroke-width="0.3"/>'
                    )
            parts.append(
                f'<text x="{origin_x + panel_width / 2}" y="{origin_y + k * cell + 18}" '
                f'text-anchor="middle" font-family="sans-serif" font-size="11">Mesh X</text>'
            )
            parts.append(
                f'<text x="{origin_x - 24}" y="{origin_y + k * cell / 2}" '
                f'text-anchor="middle" font-family="sans-serif" font-size="11" '
                f'transform="rotate(-90 {origin_x - 24} {origin_y + k * cell / 2})">Mesh Y</text>'
            )
        parts.append("</svg>")
        (output / f"heatmap_r{ranks}_{profile}_fold_{fold}.svg").write_text(
            "\n".join(parts) + "\n", encoding="utf-8"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--output", type=Path, default=ROOT / "experiments/wse_phase3"
    )
    parser.add_argument("--compute-capability", type=float, default=2.5e10)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    results = []
    all_stages = []
    all_links = []
    for ranks, case_dir in DEFAULT_CASES.items():
        plan = case_dir / "wse_plan.json"
        cost = case_dir / "static_cost_estimate.json"
        for profile in ("wse_fast", "iq"):
            for fold in ("on", "off"):
                tag = f"r{ranks}_{profile}_fold_{fold}"
                prefix = output / tag
                command = [
                    sys.executable,
                    str(ROOT / "wse_compiler.py"),
                    str(plan),
                    str(cost),
                    str(ROOT / "booksim2/ccdg_lammps_4x4.cfg"),
                    "--compute-capability",
                    str(args.compute_capability),
                    "-o",
                    str(prefix),
                ]
                if profile == "wse_fast":
                    command.append("--wse-fast-profile")
                if fold == "off":
                    command.append("--no-fold")
                run(command, output / f"{tag}.compile.log")
                run(
                    [
                        str(ROOT / "booksim2/run_wse_program.sh"),
                        str(prefix.with_suffix(".program.json")),
                        str(ROOT / "booksim2/ccdg_lammps_4x4.cfg"),
                        str(args.timeout),
                    ],
                    output / f"{tag}.acceptance.log",
                )
                program = read_json(prefix.with_suffix(".program.json"))
                report = read_json(prefix.with_suffix(".report.json"))
                acceptance = read_json(prefix.with_suffix(".acceptance.json"))
                replay = read_json(prefix.with_suffix(".replay.ccdg"))
                stats_text = prefix.with_suffix(".booksim.stats").read_text(encoding="utf-8")
                result = {
                    "ranks": ranks,
                    "profile": profile,
                    "fold": fold,
                    "compiler_cycles": acceptance["compiler_cycles"],
                    "booksim_cycles": acceptance["booksim_cycles"],
                    "compute_barrier_cycles": report["cycles"]["compute_barrier"],
                    "steps_per_sec": (
                        program["hardware"]["noc_frequency_ghz"]
                        * 1e9
                        / acceptance["booksim_cycles"]
                    ),
                    "relative_error": acceptance["relative_error"],
                    "congestion_ratio": acceptance["congestion_ratio"],
                    "wavefronts": acceptance["wavefronts"]["injected"],
                    "commands": acceptance["commands"]["delivered"],
                    "branches": acceptance["branches"]["delivered"],
                    "logical_bytes": report["conservation"]["input_bytes"],
                    "branch_link_bytes": report["replay"]["branch_link_bytes"],
                    "status": acceptance["status"],
                }
                results.append(result)
                for row in flatten_stages(
                    program,
                    replay,
                    prefix.with_suffix(".replay.est"),
                    stats_text,
                ):
                    all_stages.append(
                        {"ranks": ranks, "profile": profile, "fold": fold, **row}
                    )
                all_links.extend(link_rows(program, ranks, profile, fold))

    baselines = historical_baselines()
    wse_lookup = {
        row["ranks"]: row
        for row in results
        if row["profile"] == "wse_fast" and row["fold"] == "on"
    }
    comparisons = []
    for baseline in baselines:
        wse = wse_lookup.get(baseline["ranks"])
        if wse:
            comparisons.append(
                {
                    **baseline,
                    "wse_cycles": wse["booksim_cycles"],
                    "wse_vs_baseline_gain_pct": (
                        (baseline["cycles"] - wse["booksim_cycles"])
                        * 100.0
                        / baseline["cycles"]
                    ),
                }
            )

    payload = {
        "schema_version": 1,
        "case": "LiAlOCl 2688 atoms, one short-range steady timestep",
        "compute_capability_ops_s": args.compute_capability,
        "results": results,
        "per_stage": all_stages,
        "baselines": baselines,
        "comparisons": comparisons,
    }
    (output / "summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    write_csv(output / "matrix.csv", results)
    write_csv(output / "per_stage.csv", all_stages)
    write_csv(output / "link_utilization.csv", all_links)
    write_csv(output / "baseline_comparison.csv", comparisons)
    render_heatmaps(output, results, all_links)
    render_svg_heatmaps(output, results, all_links)
    print(f"WSE Phase-3 matrix PASS: {len(results)} points; output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
