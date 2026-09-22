#!/usr/bin/env python3
"""Build strong-scaling tables from uncalibrated cost bounds and BookSim."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from estimate_lammps_cost import estimate, parse_system

ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Scenario:
    name: str
    input_file: Path
    plan_file: Path
    profile_file: Path


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def run(command: list[str], cwd: Path) -> str:
    result = subprocess.run(
        command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"{result.stdout}"
        )
    return result.stdout


def parse_scenario(values: list[str]) -> Scenario:
    name, input_file, plan_file, profile_file = values
    scenario = Scenario(
        name=name,
        input_file=Path(input_file).resolve(),
        plan_file=Path(plan_file).resolve(),
        profile_file=Path(profile_file).resolve(),
    )
    for path in (scenario.input_file, scenario.plan_file, scenario.profile_file):
        if not path.is_file():
            raise ValueError(f"missing scenario input: {path}")
    return scenario


def measured_compute(profile: dict[str, Any]) -> dict[str, float | None]:
    ranks = profile["ranks"]
    def metric(row: dict[str, Any], current: str, legacy: str) -> float:
        return float(row[current] if current in row else row[legacy])

    instruction_values = [
        metric(row, "retired_instructions_per_step", "instructions_per_step")
        for row in ranks
    ]
    cycle_values = [
        metric(row, "hardware_cycles_per_step", "cycles_per_step")
        for row in ranks
    ]
    wall_seconds = [
        float(row["wall_seconds_per_step"])
        for row in ranks
        if "wall_seconds_per_step" in row
    ]
    wall_cycles_nominal = [
        float(row["wall_cycles_nominal_per_step"])
        for row in ranks
        if "wall_cycles_nominal_per_step" in row
    ]
    dp_ops = profile.get("dp_ops_per_step", {}).get("sum")
    return {
        "measured_retired_instructions_per_step": sum(instruction_values),
        "measured_retired_instructions_per_rank_average": (
            sum(instruction_values) / len(instruction_values)
        ),
        "measured_retired_instructions_per_rank_maximum": max(instruction_values),
        "measured_hardware_cycles_per_step_sum": sum(cycle_values),
        "measured_hardware_cycles_per_rank_maximum": max(cycle_values),
        "measured_critical_wall_seconds_per_step": max(wall_seconds) if wall_seconds else None,
        "measured_critical_wall_cycles_nominal_per_step": (
            max(wall_cycles_nominal) if wall_cycles_nominal else None
        ),
        "measured_dp_ops_per_step_optional": float(dp_ops) if dp_ops is not None else None,
    }


def analyze_one(
    scenario: Scenario,
    cfg: Path,
    compute_capability: float,
    timeout: int,
    workdir: Path,
) -> dict[str, Any]:
    static = estimate(parse_system(scenario.input_file))
    static["input_file"] = str(scenario.input_file)
    ranks = int(static["input"]["num_ranks"])
    profile = read_json(scenario.profile_file)
    if int(profile["num_ranks"]) != ranks:
        raise ValueError(
            f"{scenario.name}: profile ranks={profile['num_ranks']} != input ranks={ranks}"
        )
    profile_input = profile.get("input_file")
    if profile_input and Path(profile_input).resolve() != scenario.input_file:
        raise ValueError(f"{scenario.name}: compute profile belongs to another input")

    prefix = workdir / scenario.name
    static_path = prefix.with_suffix(".static.json")
    static_path.write_text(json.dumps(static, indent=2) + "\n", encoding="utf-8")
    compile_log = run(
        [
            "python3",
            str(ROOT / "wse_compiler.py"),
            str(scenario.plan_file),
            str(static_path),
            str(cfg),
            "--wse-fast-profile",
            "--compute-capability",
            str(compute_capability),
            "-o",
            str(prefix),
        ],
        ROOT,
    )
    replay_log = run(
        [
            str(ROOT / "booksim2/run_wse_program.sh"),
            str(prefix.with_suffix(".program.json")),
            str(cfg),
            str(timeout),
        ],
        ROOT,
    )
    report = read_json(prefix.with_suffix(".report.json"))
    acceptance = read_json(prefix.with_suffix(".acceptance.json"))
    expected_total = int(acceptance["compiler_cycles"])
    expected_compute = int(report["cycles"]["compute_barrier"])
    expected_communication = expected_total - expected_compute
    compute_to_communication = (
        expected_compute / expected_communication
        if expected_communication > 0
        else math.inf
    )
    communication_to_compute = (
        expected_communication / expected_compute if expected_compute > 0 else math.inf
    )
    summary = static["summary"]
    mesh = int(round(math.sqrt(ranks)))
    row: dict[str, Any] = {
        "scenario": scenario.name,
        "atoms": int(static["input"]["atoms"]),
        "ranks": ranks,
        "mesh": f"{mesh}x{mesh}",
        "scenario_scale": f"{int(static['input']['atoms'])} atoms / {ranks} ranks ({mesh}x{mesh})",
        "compute_to_communication_ratio": compute_to_communication,
        "communication_to_compute_ratio": communication_to_compute,
        "expected_compute_cycles": expected_compute,
        "expected_communication_cycles": expected_communication,
        "expected_cycles": expected_total,
        "booksim_measured_cycles": int(acceptance["booksim_cycles"]),
        "booksim_relative_error": float(acceptance["relative_error"]),
        "theoretical_compute_ops_lower": float(
            summary["C1_steady_ops_lower"]["sum"]
        ),
        "theoretical_compute_ops_upper": float(
            summary["C1_steady_ops_upper"]["sum"]
        ),
        "theoretical_critical_rank_ops_lower": float(
            summary["C1_steady_ops_lower"]["maximum"]
        ),
        "theoretical_critical_rank_ops_upper": float(
            summary["C1_steady_ops_upper"]["maximum"]
        ),
        "theoretical_communication_bytes_lower": float(
            summary["T1_send_bytes_lower"]["sum"]
        ),
        "theoretical_communication_bytes_upper": float(
            summary["T1_send_bytes_upper"]["sum"]
        ),
        **measured_compute(profile),
        "compute_profile_steps": int(profile["steps"]),
        "compute_profile_repeats": int(profile["repeats"]),
        "compute_profile_calibration": profile.get("calibration", "legacy-none"),
        "compute_profile_host_logical_cpus": profile.get(
            "host_logical_cpus", os.cpu_count()
        ),
        "compute_profile_oversubscribed": profile.get(
            "oversubscribed", bool(os.cpu_count() and ranks > os.cpu_count())
        ),
        "compute_profile_cpu_frequency": profile.get("cpu_frequency"),
        "booksim_average_packet_queue_cycles": acceptance.get(
            "average_packet_queue_cycles"
        ),
        "booksim_average_injection_rate": acceptance.get("average_injection_rate"),
        "booksim_injection_saturation_ratio": acceptance.get(
            "injection_saturation_ratio"
        ),
        "booksim_saturated_injection_rate": acceptance.get(
            "saturated_injection_rate"
        ),
        "booksim_aggregate_communication_to_compute_ratio": acceptance.get(
            "communication_to_compute_ratio"
        ),
        "input_file": str(scenario.input_file),
        "plan_file": str(scenario.plan_file),
        "profile_method": profile["method"],
        "compile_status": compile_log.strip().splitlines()[0],
        "replay_status": replay_log.strip().splitlines()[-1],
    }
    return row


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def render_markdown(rows: list[dict[str, Any]]) -> str:
    lines = [
        "# 强扩展计算通信分析",
        "",
        "口径：WSE-fast、fold on；计算通信比 = 预期计算 cycles / 预期通信 cycles；",
        "预期 cycles = 编译器计算 barrier + 编译通信窗口。未使用 profile 校准系数。",
        "",
        "| 场景规模 | 计算通信比 | 预期 cycles |",
        "|---|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['scenario_scale']} | "
            f"{row['compute_to_communication_ratio']:.4f} | "
            f"{row['expected_cycles']:,} |"
        )
    lines.extend([
        "",
        "## 估算与实测计算量",
        "",
        "理论量是源码推导的 algorithmic ops 区间；实测量是 retired instructions。",
        "二者单位不同，仅并列观察趋势，不用经验系数互相换算。",
        "",
        "| Ranks | 理论总 C1 ops 区间 | 实测 DP ops/step | "
        "实测总 instructions/step | 预期 compute/comm cycles | BookSim cycles | 误差 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in rows:
        lines.append(
            f"| {row['ranks']} | "
            f"{row['theoretical_compute_ops_lower']:,.0f}–"
            f"{row['theoretical_compute_ops_upper']:,.0f} | "
            f"{row['measured_dp_ops_per_step_optional']:,.0f} | "
            f"{row['measured_retired_instructions_per_step']:,.0f} | "
            f"{row['expected_compute_cycles']:,}/"
            f"{row['expected_communication_cycles']:,} | "
            f"{row['booksim_measured_cycles']:,} | "
            f"{row['booksim_relative_error']:.4%} |"
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scenario",
        action="append",
        nargs=4,
        metavar=("NAME", "INPUT", "PLAN", "PROFILE"),
        required=True,
    )
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=ROOT / "experiments/scaling"
    )
    parser.add_argument(
        "--booksim-config",
        type=Path,
        default=ROOT / "booksim2/ccdg_lammps_4x4.cfg",
    )
    parser.add_argument("--compute-capability", type=float, default=2.5e10)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    scenarios = [parse_scenario(values) for values in args.scenario]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="dpmd-scaling-") as tmp:
        rows = [
            analyze_one(
                scenario,
                args.booksim_config.resolve(),
                args.compute_capability,
                args.timeout,
                Path(tmp),
            )
            for scenario in scenarios
        ]
    rows.sort(key=lambda row: int(row["ranks"]))
    base = rows[0]
    for row in rows:
        row["theoretical_critical_ops_relative_to_smallest_rank"] = (
            row["theoretical_critical_rank_ops_upper"]
            / base["theoretical_critical_rank_ops_upper"]
        )
        row["measured_total_instructions_relative_to_smallest_rank"] = (
            row["measured_retired_instructions_per_step"]
            / base["measured_retired_instructions_per_step"]
        )
        row["measured_max_rank_instructions_relative_to_smallest_rank"] = (
            row["measured_retired_instructions_per_rank_maximum"]
            / base["measured_retired_instructions_per_rank_maximum"]
        )
    summary_fields = [
        "scenario_scale",
        "compute_to_communication_ratio",
        "expected_cycles",
    ]
    write_csv(output / "scaling_summary.csv", rows, summary_fields)
    write_csv(output / "scaling_detail.csv", rows, list(rows[0]))
    (output / "scaling_summary.md").write_text(
        render_markdown(rows), encoding="utf-8"
    )
    metadata = {
        "schema_version": 1,
        "calibration": "none",
        "compute_capability_ops_s": args.compute_capability,
        "profile": "wse_fast",
        "fold_pbc": True,
        "definitions": {
            "compute_to_communication_ratio": (
                "critical-path expected_compute_cycles / expected_communication_cycles"
            ),
            "expected_cycles": (
                "WSE compiler manager-replay total using theoretical C1 upper bound"
            ),
            "calculation_comparison": (
                "theoretical algorithmic ops and measured retired instructions are "
                "reported in separate units without fitted conversion"
            ),
            "booksim_aggregate_communication_to_compute_ratio": (
                "sum over ranks of blocked+congestion+sched_wait cycles divided "
                "by sum over ranks of compute dwell; not the inverse of the "
                "critical-path ratio"
            ),
        },
        "rows": rows,
    }
    (output / "scaling_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(render_markdown(rows), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
