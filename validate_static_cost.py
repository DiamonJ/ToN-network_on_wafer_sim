#!/usr/bin/env python3
"""Compare uncalibrated static T1/C1 bounds with captured evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def captured_t1(plan: dict[str, Any]) -> tuple[list[float], list[float]]:
    ranks = int(plan["num_ranks"])
    steady = [0.0] * ranks
    rebuild_extra = [0.0] * ranks
    rebuild_phases = {"borders", "reverse", "pair_forward", "pair_reverse"}
    for record in plan["records"]:
        if record.get("kind") != "message":
            continue
        rank, size = int(record["rank"]), float(record["bytes"])
        if record.get("scope") == "run":
            steady[rank] += size
        elif (
            record.get("scope") == "setup"
            and record.get("phase") in rebuild_phases
        ):
            rebuild_extra[rank] += size
    return steady, rebuild_extra


def interval_report(
    lower: list[float], upper: list[float], observed: list[float]
) -> dict[str, Any]:
    if not lower or len(lower) != len(upper) or len(lower) != len(observed):
        raise ValueError("bound and observed rank arrays must be non-empty and equal length")
    if any(lo < 0 or hi < lo for lo, hi in zip(lower, upper)):
        raise ValueError("every theoretical interval must satisfy 0 <= lower <= upper")
    gaps = []
    inside = []
    for lo, hi, value in zip(lower, upper, observed):
        ok = lo <= value <= hi
        inside.append(ok)
        gaps.append(0.0 if ok else lo - value if value < lo else value - hi)
    lower_sum, upper_sum, observed_sum = sum(lower), sum(upper), sum(observed)
    return {
        "theoretical_lower_sum": lower_sum,
        "theoretical_upper_sum": upper_sum,
        "observed_sum": observed_sum,
        "global_inside_interval": lower_sum <= observed_sum <= upper_sum,
        "ranks_inside_interval": sum(inside),
        "rank_count": len(observed),
        "all_ranks_inside_interval": all(inside),
        "per_rank_inside_interval": inside,
        "per_rank_distance_to_interval": gaps,
        "maximum_distance_to_interval": max(gaps),
    }


def distribution(values: list[float]) -> dict[str, float]:
    return {
        "sum": sum(values),
        "average": sum(values) / len(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def profile_field(row: dict[str, Any], new: str, legacy: str) -> float:
    if new in row:
        return float(row[new])
    if legacy in row:
        return float(row[legacy])
    raise ValueError(f"compute profile lacks {new!r}")


def compare(
    estimate: dict[str, Any],
    plan: dict[str, Any],
    profile: dict[str, Any],
) -> dict[str, Any]:
    ranks = int(estimate["input"]["num_ranks"])
    if int(plan["num_ranks"]) != ranks or int(profile["num_ranks"]) != ranks:
        raise ValueError("rank count differs between estimate, plan, and profile")
    estimated = sorted(estimate["ranks"], key=lambda row: int(row["rank"]))
    measured = sorted(profile["ranks"], key=lambda row: int(row["rank"]))
    expected_ranks = list(range(ranks))
    if [int(row["rank"]) for row in estimated] != expected_ranks:
        raise ValueError("static estimate does not contain every rank exactly once")
    if [int(row["rank"]) for row in measured] != expected_ranks:
        raise ValueError("compute profile does not contain every rank exactly once")
    estimate_input = estimate.get("input_file")
    profile_input = profile.get("input_file")
    if estimate_input and profile_input:
        if Path(estimate_input).resolve() != Path(profile_input).resolve():
            raise ValueError(
                "static estimate and compute profile refer to different input files"
            )

    observed_steady, observed_rebuild_extra = captured_t1(plan)
    theoretical_steady_lower = [
        float(row.get("T1_steady_send_bytes_lower", row["T1_steady_send_bytes"]))
        for row in estimated
    ]
    theoretical_steady_upper = [
        float(row.get("T1_steady_send_bytes_upper", row["T1_steady_send_bytes"]))
        for row in estimated
    ]
    theoretical_rebuild_lower = [
        float(row.get("T1_rebuild_send_bytes_lower", row["T1_rebuild_send_bytes"]))
        for row in estimated
    ]
    theoretical_rebuild_upper = [
        float(row.get("T1_rebuild_send_bytes_upper", row["T1_rebuild_send_bytes"]))
        for row in estimated
    ]
    theoretical_lower = [
        min(steady, rebuild)
        for steady, rebuild in zip(theoretical_steady_lower, theoretical_rebuild_lower)
    ]
    theoretical_upper = [
        max(steady, rebuild)
        for steady, rebuild in zip(theoretical_steady_upper, theoretical_rebuild_upper)
    ]
    t1_steady = interval_report(
        theoretical_steady_lower, theoretical_steady_upper, observed_steady
    )
    t1_rebuild = interval_report(
        theoretical_rebuild_lower, theoretical_rebuild_upper, observed_rebuild_extra
    )
    t1 = {
        "theoretical_lower_sum": sum(theoretical_lower),
        "theoretical_upper_sum": sum(theoretical_upper),
        "observed_lower_sum": min(sum(observed_steady), sum(observed_rebuild_extra)),
        "observed_upper_sum": max(sum(observed_steady), sum(observed_rebuild_extra)),
        "steady": t1_steady,
        "rebuild_without_exchange": t1_rebuild,
        "all_observed_scenarios_inside_interval": (
            t1_steady["global_inside_interval"]
            and t1_rebuild["global_inside_interval"]
        ),
        "semantics": (
            "bounds span a normal step and a neighbor-rebuild step; these "
            "alternative scenarios are not added together"
        ),
    }

    c1_lower = [
        float(row.get("C1_steady_ops_lower", row["C1_steady_ops_min"]))
        for row in estimated
    ]
    c1_upper = [
        float(row.get("C1_steady_ops_upper", row["C1_steady_ops_max"]))
        for row in estimated
    ]
    instructions = [
        profile_field(row, "retired_instructions_per_step", "instructions_per_step")
        for row in measured
    ]
    hardware_cycles = [
        profile_field(row, "hardware_cycles_per_step", "cycles_per_step")
        for row in measured
    ]
    wall_seconds = [
        float(row["wall_seconds_per_step"])
        for row in measured
        if "wall_seconds_per_step" in row
    ]
    measured_c1: dict[str, Any] = {
        "retired_instructions_per_step": distribution(instructions),
        "hardware_cycles_per_step": distribution(hardware_cycles),
        "hardware_cpi_per_rank": [
            cycles / insn for cycles, insn in zip(hardware_cycles, instructions)
        ],
        "static_algorithmic_ops_lower": distribution(c1_lower),
        "static_algorithmic_ops_upper": distribution(c1_upper),
        "comparison_policy": (
            "reported side by side only: algorithmic operations and retired "
            "instructions are different units; no fitted conversion coefficient"
        ),
    }
    if wall_seconds:
        if len(wall_seconds) != ranks:
            raise ValueError("wall time must be available for all ranks or none")
        wall_lower = [float(row["wall_cycles_lower_per_step"]) for row in measured]
        wall_nominal = [float(row["wall_cycles_nominal_per_step"]) for row in measured]
        wall_upper = [float(row["wall_cycles_upper_per_step"]) for row in measured]
        measured_c1.update({
            "wall_seconds_per_step": distribution(wall_seconds),
            "wall_cycles_lower_per_step": distribution(wall_lower),
            "wall_cycles_nominal_per_step": distribution(wall_nominal),
            "wall_cycles_upper_per_step": distribution(wall_upper),
            "hardware_cycles_inside_wall_frequency_bounds": all(
                lo <= cycles <= hi
                for lo, cycles, hi in zip(wall_lower, hardware_cycles, wall_upper)
            ),
            "retired_instructions_per_wall_cycle_nominal_per_rank": [
                insn / cycles
                for insn, cycles in zip(instructions, wall_nominal)
            ],
            "retired_instructions_per_wall_cycle_bounds_per_rank": [
                [insn / hi, insn / lo]
                for insn, lo, hi in zip(instructions, wall_lower, wall_upper)
            ],
        })

    static_bounds_valid = all(lo <= hi for lo, hi in zip(c1_lower, c1_upper))
    evidence_complete = all(value > 0 for value in instructions + hardware_cycles)
    passed = (
        t1["all_observed_scenarios_inside_interval"]
        and static_bounds_valid
        and evidence_complete
    )
    return {
        "schema_version": 2,
        "validation_policy": "uncalibrated_bounds",
        "calibration": "none",
        "num_ranks": ranks,
        "input_file": estimate.get("input_file"),
        "input": estimate["input"],
        "profile_method": profile["method"],
        "profile_event_encoding": profile["event_encoding"],
        "profile_steps": profile["steps"],
        "cpu_frequency": profile.get("cpu_frequency"),
        "t1_iteration_send_bytes": t1,
        "c1": measured_c1,
        "passed": passed,
        "gates": {
            "t1_scenarios_inside_theoretical_bounds": (
                t1["all_observed_scenarios_inside_interval"]
            ),
            "c1_static_bounds_well_formed": static_bounds_valid,
            "c1_measurement_complete": evidence_complete,
        },
    }


def render(result: dict[str, Any]) -> str:
    data, t1, c1 = result["input"], result["t1_iteration_send_bytes"], result["c1"]
    lines = [
        f"{data['atoms']} atoms / {result['num_ranks']} ranks",
        "",
        "T1 理论通信量边界（send payload bytes/step）：",
        f"理论区间：[{t1['theoretical_lower_sum']:,.0f}, {t1['theoretical_upper_sum']:,.0f}]",
        (
            f"源码插桩区间：[{t1['observed_lower_sum']:,.0f}, "
            f"{t1['observed_upper_sum']:,.0f}]"
        ),
        (
            "steady/rebuild 均落入区间："
            f"{'YES' if t1['all_observed_scenarios_inside_interval'] else 'NO'}"
        ),
        "",
        "C1 动态测量（不做系数校准）：",
        (
            "退休指令："
            f"{c1['retired_instructions_per_step']['sum']:,.0f} instructions/step"
        ),
        f"硬件周期：{c1['hardware_cycles_per_step']['sum']:,.0f} cycles/step",
    ]
    if "wall_seconds_per_step" in c1:
        lines.extend([
            f"墙上时间：{c1['wall_seconds_per_step']['maximum']:.9f} s/step (critical rank)",
            (
                "墙时×CPU频率周期区间："
                f"[{c1['wall_cycles_lower_per_step']['maximum']:,.0f}, "
                f"{c1['wall_cycles_upper_per_step']['maximum']:,.0f}] cycles/step"
            ),
        ])
    lines.extend([
        "",
        "C1 静态理论算法操作数边界（与退休指令不同单位）：",
        (
            f"[{c1['static_algorithmic_ops_lower']['sum']:,.0f}, "
            f"{c1['static_algorithmic_ops_upper']['sum']:,.0f}] ops/step"
        ),
        "未使用 profile 反推或拟合任何换算系数。",
        "",
        f"OVERALL: {'PASS' if result['passed'] else 'FAIL'}",
    ])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("estimate", type=Path)
    parser.add_argument("wse_plan", type=Path)
    parser.add_argument("compute_profile", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args()
    result = compare(
        read_json(args.estimate),
        read_json(args.wse_plan),
        read_json(args.compute_profile),
    )
    text = render(result)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.with_suffix(".json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
        args.output.with_suffix(".txt").write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
