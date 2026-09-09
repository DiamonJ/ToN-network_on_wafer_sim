#!/usr/bin/env python3
"""Compare source-level CommBrick run traffic with a trim-only CCDG."""

from __future__ import annotations

import argparse
import collections
import json
import math
from typing import Any


Direction = tuple[int, int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("plan")
    parser.add_argument("ccdg")
    parser.add_argument("--threshold", type=float, default=0.02)
    return parser.parse_args()


def direction(src: int, dst: int, k: int) -> Direction:
    return (dst % k - src % k, dst // k - src // k)


def format_direction(value: Direction) -> str:
    return f"{value[0]},{value[1]}"


def main() -> int:
    args = parse_args()
    with open(args.plan, encoding="utf-8") as stream:
        plan: dict[str, Any] = json.load(stream)
    with open(args.ccdg, encoding="utf-8") as stream:
        ccdg: dict[str, Any] = json.load(stream)

    ranks = int(plan["num_ranks"])
    if int(ccdg["num_ranks"]) != ranks:
        raise SystemExit(
            f"ERROR: rank mismatch: plan={ranks} ccdg={ccdg['num_ranks']}"
        )
    k = int(round(math.sqrt(ranks)))
    if k * k != ranks:
        raise SystemExit(f"ERROR: num_ranks={ranks} is not a square mesh")

    plan_bytes: collections.Counter[Direction] = collections.Counter()
    plan_messages: collections.Counter[Direction] = collections.Counter()
    phase_bytes: collections.Counter[str] = collections.Counter()
    plan_collectives: collections.Counter[tuple[str, int]] = collections.Counter()
    for record in plan.get("records", []):
        if record.get("scope") != "run":
            continue
        if record.get("kind") == "collective":
            key = (str(record.get("operation", "")).upper(), int(record["bytes"]))
            plan_collectives[key] += 1
            continue
        if record.get("kind") != "message":
            continue
        src = int(record["src"])
        dst = int(record["dst"])
        key = direction(src, dst, k)
        size = int(record.get("bytes", 0))
        plan_bytes[key] += size
        plan_messages[key] += 1
        phase_bytes[str(record.get("phase", "unknown"))] += size

    ccdg_bytes: collections.Counter[Direction] = collections.Counter()
    ccdg_messages: collections.Counter[Direction] = collections.Counter()
    ccdg_collectives: collections.Counter[tuple[str, int]] = collections.Counter()
    for node in ccdg.get("nodes", []):
        if node.get("type") in ("ALLREDUCE", "ALLTOALLV"):
            ccdg_collectives[(str(node["type"]), int(node.get("comm_bytes", 0)))] += 1
        if node.get("type") not in ("SEND", "ISEND"):
            continue
        src = int(node["rank"])
        dst = int(node.get("comm_dst", src))
        key = direction(src, dst, k)
        ccdg_bytes[key] += int(node.get("comm_bytes", 0))
        ccdg_messages[key] += 1

    comparisons = []
    passed = True
    for key in sorted(set(plan_bytes) | set(ccdg_bytes)):
        source_bytes = int(plan_bytes[key])
        trace_bytes = int(ccdg_bytes[key])
        if trace_bytes:
            relative_error = abs(source_bytes - trace_bytes) / trace_bytes
        else:
            relative_error = 0.0 if source_bytes == 0 else float("inf")
        direction_passed = relative_error <= args.threshold
        passed = passed and direction_passed
        comparisons.append(
            {
                "direction": format_direction(key),
                "plan_bytes": source_bytes,
                "ccdg_bytes": trace_bytes,
                "plan_messages": int(plan_messages[key]),
                "ccdg_messages": int(ccdg_messages[key]),
                "relative_error": relative_error,
                "passed": direction_passed,
            }
        )

    collective_comparisons = []
    collectives_passed = True
    for key in sorted(plan_collectives):
        expected_count = int(plan_collectives[key])
        observed_count = int(ccdg_collectives[key])
        item_passed = observed_count >= expected_count
        collectives_passed = collectives_passed and item_passed
        collective_comparisons.append(
            {
                "operation": key[0],
                "bytes": key[1],
                "plan_count": expected_count,
                "ccdg_count": observed_count,
                "passed": item_passed,
            }
        )

    report = {
        "passed": passed and bool(comparisons) and collectives_passed,
        "threshold": args.threshold,
        "num_ranks": ranks,
        "plan_total_bytes": sum(plan_bytes.values()),
        "ccdg_total_bytes": sum(ccdg_bytes.values()),
        "plan_total_messages": sum(plan_messages.values()),
        "ccdg_total_messages": sum(ccdg_messages.values()),
        "plan_phase_bytes": dict(sorted(phase_bytes.items())),
        "directions": comparisons,
        "kspace_collectives_covered": collectives_passed,
        "kspace_collectives": collective_comparisons,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
