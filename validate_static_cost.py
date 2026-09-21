#!/usr/bin/env python3
"""Compare static T1/C1 estimates with LAMMPS communication and PMU profiles."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def stats(estimate: list[float], actual: list[float]) -> dict[str, Any]:
    if len(estimate) != len(actual) or not estimate:
        raise ValueError("estimate/actual rank arrays must be non-empty and equal length")
    if any(value <= 0 for value in actual):
        raise ValueError("actual per-rank values must be positive")
    est_sum, actual_sum = sum(estimate), sum(actual)
    errors = [(est - obs) / obs for est, obs in zip(estimate, actual)]
    return {
        "estimate_sum": est_sum,
        "actual_sum": actual_sum,
        "global_relative_error": (est_sum - actual_sum) / actual_sum,
        "max_rank_absolute_relative_error": max(abs(value) for value in errors),
        "mean_rank_absolute_relative_error": sum(abs(value) for value in errors)
        / len(errors),
        "per_rank_relative_error": errors,
    }


def captured_t1(plan: dict[str, Any]) -> tuple[list[float], list[float]]:
    ranks = int(plan["num_ranks"])
    steady = [0.0] * ranks
    rebuild = [0.0] * ranks
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
            rebuild[rank] += size
    return steady, rebuild


def compare(
    estimate: dict[str, Any],
    plan: dict[str, Any],
    profile: dict[str, Any],
    t1_threshold: float,
    c1_global_threshold: float,
    c1_rank_threshold: float,
) -> dict[str, Any]:
    ranks = int(estimate["input"]["num_ranks"])
    if int(plan["num_ranks"]) != ranks or int(profile["num_ranks"]) != ranks:
        raise ValueError("rank count differs between estimate, plan, and profile")
    estimated_ranks = sorted(estimate["ranks"], key=lambda row: int(row["rank"]))
    profiled_ranks = sorted(profile["ranks"], key=lambda row: int(row["rank"]))
    if [int(row["rank"]) for row in estimated_ranks] != list(range(ranks)):
        raise ValueError("static estimate does not contain every rank exactly once")
    if [int(row["rank"]) for row in profiled_ranks] != list(range(ranks)):
        raise ValueError("PMU profile does not contain every rank exactly once")

    actual_steady, actual_rebuild = captured_t1(plan)
    t1_steady = stats(
        [float(row["T1_steady_send_bytes"]) for row in estimated_ranks],
        actual_steady,
    )
    t1_rebuild = stats(
        [float(row["T1_rebuild_send_bytes"]) for row in estimated_ranks],
        actual_rebuild,
    )
    c1 = stats(
        [float(row["C1_steady_ops_midpoint"]) for row in estimated_ranks],
        [float(row["dp_ops_per_step"]) for row in profiled_ranks],
    )
    c1_min = [float(row["C1_steady_ops_min"]) for row in estimated_ranks]
    c1_max = [float(row["C1_steady_ops_max"]) for row in estimated_ranks]
    c1_actual = [float(row["dp_ops_per_step"]) for row in profiled_ranks]
    c1["estimate_min_sum"] = sum(c1_min)
    c1["estimate_max_sum"] = sum(c1_max)
    c1["ranks_inside_estimate_interval"] = sum(
        low <= observed <= high
        for low, high, observed in zip(c1_min, c1_max, c1_actual)
    )
    c1["rank_count"] = ranks

    t1_passed = (
        abs(t1_steady["global_relative_error"]) <= t1_threshold
        and t1_steady["max_rank_absolute_relative_error"] <= t1_threshold
        and abs(t1_rebuild["global_relative_error"]) <= t1_threshold
        and t1_rebuild["max_rank_absolute_relative_error"] <= t1_threshold
    )
    c1_passed = (
        abs(c1["global_relative_error"]) <= c1_global_threshold
        and c1["max_rank_absolute_relative_error"] <= c1_rank_threshold
    )
    return {
        "schema_version": 1,
        "num_ranks": ranks,
        "input_file": estimate.get("input_file"),
        "input": estimate["input"],
        "profile_method": profile["method"],
        "profile_event_encoding": profile["event_encoding"],
        "profile_steps": profile["steps"],
        "thresholds": {
            "t1_absolute_relative_error": t1_threshold,
            "c1_global_absolute_relative_error": c1_global_threshold,
            "c1_rank_absolute_relative_error": c1_rank_threshold,
        },
        "t1_steady": t1_steady,
        "t1_rebuild_without_exchange": t1_rebuild,
        "c1_steady_dp_ops": c1,
        "passed": t1_passed and c1_passed,
        "gates": {"t1": t1_passed, "c1": c1_passed},
    }


def percent(value: float) -> str:
    return f"{value * 100:.3f}%"


def render(result: dict[str, Any]) -> str:
    data = result["input"]
    t1 = result["t1_steady"]
    rebuild = result["t1_rebuild_without_exchange"]
    c1 = result["c1_steady_dp_ops"]
    input_path = str(result.get("input_file") or "").lower()
    system = (
        "LiAlOCl"
        if "lialocl" in input_path
        else "Cu"
        if "_cu_" in input_path
        else "H2O"
        if "h2o" in input_path
        else data["pair_style"]
    )
    lines = [
        (
            f"{system} {data['atoms']} atoms / "
            f"{result['num_ranks']} ranks"
        ),
        "",
        "稳态 T1：",
        f"估算器：{t1['estimate_sum']:,.0f} bytes/step",
        f"真实捕捉：{t1['actual_sum']:,.0f} bytes/step",
        f"全局误差：{percent(t1['global_relative_error'])}",
        f"逐 rank 最大绝对误差：{percent(t1['max_rank_absolute_relative_error'])}",
        "",
        "重建 T1（borders + reverse，不含 exchange）：",
        f"估算器：{rebuild['estimate_sum']:,.0f} bytes",
        f"真实插桩：{rebuild['actual_sum']:,.0f} bytes",
        f"全局误差：{percent(rebuild['global_relative_error'])}",
        (
            "逐 rank 最大绝对误差："
            f"{percent(rebuild['max_rank_absolute_relative_error'])}"
        ),
        "",
        "稳态 C1（double-precision arithmetic ops）：",
        f"估算器中点：{c1['estimate_sum']:,.0f} DP ops/step",
        (
            f"估算器区间：[{c1['estimate_min_sum']:,.0f}, "
            f"{c1['estimate_max_sum']:,.0f}] DP ops/step"
        ),
        f"真实 PMU：{c1['actual_sum']:,.0f} DP ops/step",
        f"全局中点误差：{percent(c1['global_relative_error'])}",
        (
            "逐 rank 最大绝对误差："
            f"{percent(c1['max_rank_absolute_relative_error'])}"
        ),
        (
            "落入估算区间的 ranks："
            f"{c1['ranks_inside_estimate_interval']}/{c1['rank_count']}"
        ),
        "",
        f"OVERALL: {'PASS' if result['passed'] else 'FAIL'}",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("estimate", type=Path)
    parser.add_argument("wse_plan", type=Path)
    parser.add_argument("flops_profile", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--t1-threshold", type=float, default=0.02)
    parser.add_argument("--c1-global-threshold", type=float, default=0.10)
    parser.add_argument("--c1-rank-threshold", type=float, default=0.15)
    args = parser.parse_args()
    result = compare(
        read_json(args.estimate),
        read_json(args.wse_plan),
        read_json(args.flops_profile),
        args.t1_threshold,
        args.c1_global_threshold,
        args.c1_rank_threshold,
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
