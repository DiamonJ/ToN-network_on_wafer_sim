#!/usr/bin/env python3
"""Summarize per-rank perf counters without fitted conversion factors."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any

FP_WEIGHTS = {
    "fp_arith_inst_retired.scalar_double": 1,
    "fp_arith_inst_retired.128b_packed_double": 2,
    "fp_arith_inst_retired.256b_packed_double": 4,
    "fp_arith_inst_retired.512b_packed_double": 8,
    "r01c7": 1,
    "r04c7": 2,
    "r10c7": 4,
    "r40c7": 8,
}
TIME_TO_SECONDS = {"ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1.0}


def read_perf_csv(path: Path, events: list[str]) -> dict[str, float]:
    values: dict[str, float] = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter=";"):
            event = row[2].strip().split(":", 1)[0] if len(row) >= 3 else ""
            if event not in events:
                continue
            value = row[0].strip().replace(",", "")
            if value.startswith("<"):
                raise RuntimeError(f"{path}: {event} = {value}")
            number = float(value)
            if event == "duration_time":
                unit = row[1].strip() if len(row) >= 2 else "ns"
                try:
                    number *= TIME_TO_SECONDS[unit]
                except KeyError as exc:
                    raise RuntimeError(
                        f"{path}: unsupported duration_time unit {unit!r}"
                    ) from exc
            values[event] = number
    missing = set(events) - values.keys()
    if missing:
        raise RuntimeError(f"{path}: missing events: {sorted(missing)}")
    return values


def distribution(values: list[float]) -> dict[str, float]:
    return {
        "sum": sum(values),
        "average": sum(values) / len(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def summarize(
    out: Path,
    ranks: int,
    steps: int,
    repeats: int,
    frequency_nominal_ghz: float,
    frequency_min_ghz: float,
    frequency_max_ghz: float,
    frequency_source: str,
    event_encoding: str,
    events: list[str],
    input_file: Path | None = None,
) -> dict[str, Any]:
    if ranks < 1 or steps < 1 or repeats < 1:
        raise ValueError("ranks, steps, and repeats must be positive")
    if not (0 < frequency_min_ghz <= frequency_nominal_ghz <= frequency_max_ghz):
        raise ValueError("CPU frequencies must satisfy 0 < min <= nominal <= max")
    required = {"cycles", "instructions", "duration_time"}
    if not required.issubset(events):
        raise ValueError(f"required perf events missing: {sorted(required - set(events))}")

    fp_events = [event for event in events if event in FP_WEIGHTS]
    records = []
    for rank in range(ranks):
        samples = []
        for repeat in range(repeats):
            base = read_perf_csv(
                out / f"baseline.repeat{repeat}.rank{rank:04d}.csv", events
            )
            run = read_perf_csv(
                out / f"run.repeat{repeat}.rank{rank:04d}.csv", events
            )
            delta = {event: run[event] - base[event] for event in events}
            if any(delta[event] <= 0 for event in required):
                raise RuntimeError(
                    f"rank {rank} repeat {repeat}: non-positive baseline-subtracted "
                    "cycles, instructions, or wall time; increase STEPS"
                )
            wall_seconds = delta["duration_time"] / steps
            sample: dict[str, Any] = {
                "repeat": repeat,
                "retired_instructions_per_step": delta["instructions"] / steps,
                "hardware_cycles_per_step": delta["cycles"] / steps,
                "wall_seconds_per_step": wall_seconds,
                "wall_cycles_nominal_per_step": (
                    wall_seconds * frequency_nominal_ghz * 1e9
                ),
                "wall_cycles_lower_per_step": wall_seconds * frequency_min_ghz * 1e9,
                "wall_cycles_upper_per_step": wall_seconds * frequency_max_ghz * 1e9,
                "baseline_events": base,
                "run_events": run,
                "delta_events": delta,
            }
            if fp_events:
                sample["dp_ops_per_step"] = sum(
                    delta[event] * FP_WEIGHTS[event] for event in fp_events
                ) / steps
            samples.append(sample)

        def median(field: str) -> float:
            return statistics.median(float(sample[field]) for sample in samples)

        row: dict[str, Any] = {
            "rank": rank,
            "steps": steps,
            "repeats": repeats,
            "aggregation": "per-rank median of baseline-subtracted repeats",
            "retired_instructions_per_step": median("retired_instructions_per_step"),
            "instructions_per_step": median("retired_instructions_per_step"),
            "hardware_cycles_per_step": median("hardware_cycles_per_step"),
            "cycles_per_step": median("hardware_cycles_per_step"),
            "wall_seconds_per_step": median("wall_seconds_per_step"),
            "wall_cycles_nominal_per_step": median("wall_cycles_nominal_per_step"),
            "wall_cycles_lower_per_step": median("wall_cycles_lower_per_step"),
            "wall_cycles_upper_per_step": median("wall_cycles_upper_per_step"),
            "repeat_records": samples,
        }
        row["C1_measured"] = {
            "retired_instructions_per_step": row["retired_instructions_per_step"],
            "hardware_cycles_per_step": row["hardware_cycles_per_step"],
            "wall_seconds_per_step": row["wall_seconds_per_step"],
            "wall_cycle_estimate_per_step": {
                "lower": row["wall_cycles_lower_per_step"],
                "nominal": row["wall_cycles_nominal_per_step"],
                "upper": row["wall_cycles_upper_per_step"],
            },
            "instructions_per_wall_cycle": {
                "lower": (
                    row["retired_instructions_per_step"]
                    / row["wall_cycles_upper_per_step"]
                ),
                "nominal": (
                    row["retired_instructions_per_step"]
                    / row["wall_cycles_nominal_per_step"]
                ),
                "upper": (
                    row["retired_instructions_per_step"]
                    / row["wall_cycles_lower_per_step"]
                ),
            },
        }
        if fp_events:
            row["dp_ops_per_step"] = median("dp_ops_per_step")
        records.append(row)

    summary: dict[str, Any] = {
        "schema_version": 2,
        "method": (
            "median of repeated perf [run(N)-run(0)]/N; retired instructions "
            "are measured directly and wall-cycle bounds equal elapsed wall time "
            "times declared CPU frequency bounds"
        ),
        "calibration": "none",
        "input_file": str(input_file.resolve()) if input_file else None,
        "num_ranks": ranks,
        "steps": steps,
        "repeats": repeats,
        "cpu_frequency": {
            "source": frequency_source,
            "nominal_ghz": frequency_nominal_ghz,
            "minimum_ghz": frequency_min_ghz,
            "maximum_ghz": frequency_max_ghz,
        },
        "event_encoding": event_encoding,
        "events": events,
        "metric_semantics": {
            "retired_instructions_per_step": "hardware instructions event",
            "hardware_cycles_per_step": "hardware unhalted core cycles event",
            "wall_seconds_per_step": "perf duration_time after run(0) subtraction",
            "wall_cycles": "wall_seconds multiplied by declared CPU frequency; no fit",
            "dp_ops_per_step": (
                "optional Intel FP_ARITH_INST_RETIRED lane-weighted diagnostic"
            ),
        },
        "ranks": records,
    }
    for field in (
        "retired_instructions_per_step",
        "hardware_cycles_per_step",
        "wall_seconds_per_step",
        "wall_cycles_nominal_per_step",
        "wall_cycles_lower_per_step",
        "wall_cycles_upper_per_step",
    ):
        summary[field] = distribution([float(row[field]) for row in records])
    if fp_events:
        summary["dp_ops_per_step"] = distribution(
            [float(row["dp_ops_per_step"]) for row in records]
        )
    summary["C1_measured"] = {
        "definition": (
            "total retired instruction work plus critical-rank wall-cycle "
            "estimate; CPU frequency is declared/measured, never fitted"
        ),
        "total_retired_instructions_per_step": summary[
            "retired_instructions_per_step"
        ]["sum"],
        "aggregate_hardware_cycles_per_step": summary[
            "hardware_cycles_per_step"
        ]["sum"],
        "critical_rank_wall_seconds_per_step": summary[
            "wall_seconds_per_step"
        ]["maximum"],
        "critical_rank_wall_cycle_estimate_per_step": {
            "lower": summary["wall_cycles_lower_per_step"]["maximum"],
            "nominal": summary["wall_cycles_nominal_per_step"]["maximum"],
            "upper": summary["wall_cycles_upper_per_step"]["maximum"],
        },
    }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile_dir", type=Path)
    parser.add_argument("ranks", type=int)
    parser.add_argument("steps", type=int)
    parser.add_argument("repeats", type=int)
    parser.add_argument("frequency_nominal_ghz", type=float)
    parser.add_argument("frequency_min_ghz", type=float)
    parser.add_argument("frequency_max_ghz", type=float)
    parser.add_argument("frequency_source")
    parser.add_argument("event_encoding")
    parser.add_argument("events", nargs="+")
    parser.add_argument("--input-file", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args()
    result = summarize(
        args.profile_dir,
        args.ranks,
        args.steps,
        args.repeats,
        args.frequency_nominal_ghz,
        args.frequency_min_ghz,
        args.frequency_max_ghz,
        args.frequency_source,
        args.event_encoding,
        args.events,
        args.input_file,
    )
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    print(json.dumps(result["retired_instructions_per_step"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
