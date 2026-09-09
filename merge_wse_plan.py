#!/usr/bin/env python3
"""Merge per-rank LAMMPS WSE-plan JSONL shards into one validated JSON file."""

from __future__ import annotations

import argparse
import collections
import glob
import json
import os
import tempfile
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", help="LAMMPS_WSE_PLAN prefix")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--expected-ranks", type=int, required=True)
    return parser.parse_args()


def read_shard(path: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    metadata: dict[str, Any] | None = None
    records: list[dict[str, Any]] = []
    with open(path, encoding="utf-8") as stream:
        for lineno, raw in enumerate(stream, 1):
            raw = raw.strip()
            if not raw:
                continue
            try:
                record = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{lineno}: invalid JSON: {exc}") from exc
            if record.get("kind") == "metadata":
                if metadata is not None:
                    raise ValueError(f"{path}: duplicate metadata record")
                metadata = record
            else:
                records.append(record)
    if metadata is None:
        raise ValueError(f"{path}: missing metadata record")
    return metadata, records


def main() -> int:
    args = parse_args()
    base_shards = sorted(glob.glob(f"{args.prefix}.rank*.jsonl"))
    component_shards = sorted(glob.glob(f"{args.prefix}.*.rank*.jsonl"))
    if len(base_shards) != args.expected_ranks:
        raise SystemExit(
            f"ERROR: expected {args.expected_ranks} WSE-plan shards, found "
            f"{len(base_shards)} for {args.prefix}.rank*.jsonl"
        )

    rank_metadata: list[dict[str, Any]] = []
    component_metadata: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []
    seen_ranks: set[int] = set()
    procgrid: tuple[int, ...] | None = None
    for shard in base_shards:
        metadata, shard_records = read_shard(shard)
        rank = int(metadata["rank"])
        if rank in seen_ranks:
            raise SystemExit(f"ERROR: duplicate rank {rank} in WSE-plan shards")
        seen_ranks.add(rank)
        if int(metadata["num_ranks"]) != args.expected_ranks:
            raise SystemExit(
                f"ERROR: rank {rank} reports num_ranks={metadata['num_ranks']}, "
                f"expected {args.expected_ranks}"
            )
        current_grid = tuple(int(v) for v in metadata["procgrid"])
        if procgrid is None:
            procgrid = current_grid
        elif current_grid != procgrid:
            raise SystemExit(
                f"ERROR: inconsistent procgrid on rank {rank}: "
                f"{current_grid} != {procgrid}"
            )
        rank_metadata.append(metadata)
        records.extend(shard_records)

    expected = set(range(args.expected_ranks))
    if seen_ranks != expected:
        raise SystemExit(
            f"ERROR: rank coverage mismatch: missing={sorted(expected-seen_ranks)} "
            f"extra={sorted(seen_ranks-expected)}"
        )

    component_ranks: dict[str, set[int]] = collections.defaultdict(set)
    for shard in component_shards:
        metadata, shard_records = read_shard(shard)
        rank = int(metadata["rank"])
        component = str(metadata.get("component", "unknown"))
        if rank in component_ranks[component]:
            raise SystemExit(
                f"ERROR: duplicate {component} shard for rank {rank}"
            )
        if int(metadata["num_ranks"]) != args.expected_ranks:
            raise SystemExit(
                f"ERROR: {component} rank {rank} reports "
                f"num_ranks={metadata['num_ranks']}, expected {args.expected_ranks}"
            )
        component_ranks[component].add(rank)
        component_metadata.append(metadata)
        records.extend(shard_records)

    for component, ranks in component_ranks.items():
        if ranks != expected:
            raise SystemExit(
                f"ERROR: {component} shard coverage mismatch: "
                f"missing={sorted(expected-ranks)} extra={sorted(ranks-expected)}"
            )

    rank_metadata.sort(key=lambda item: int(item["rank"]))
    component_metadata.sort(
        key=lambda item: (str(item.get("component", "")), int(item["rank"]))
    )
    records.sort(
        key=lambda item: (
            int(item.get("rank", -1)),
            str(item.get("component", "commbrick")),
            int(item.get("seq", -1)),
        )
    )

    summary: dict[str, dict[str, dict[str, int]]] = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: {"messages": 0, "collectives": 0, "bytes": 0}
        )
    )
    for record in records:
        kind = record.get("kind")
        if kind not in ("message", "collective"):
            continue
        scope = str(record.get("scope", "unknown"))
        phase = str(record.get("phase", "unknown"))
        if kind == "message":
            summary[scope][phase]["messages"] += 1
        else:
            summary[scope][phase]["collectives"] += 1
        summary[scope][phase]["bytes"] += int(record.get("bytes", 0))

    output = {
        "schema_version": 1,
        "num_ranks": args.expected_ranks,
        "procgrid": list(procgrid or ()),
        "rank_metadata": rank_metadata,
        "component_metadata": component_metadata,
        "summary": {scope: dict(phases) for scope, phases in summary.items()},
        "records": records,
    }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(
        prefix=f".{output_path.name}.", dir=str(output_path.parent), text=True
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(output, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temp_path, output_path)
    except BaseException:
        try:
            os.unlink(temp_path)
        except FileNotFoundError:
            pass
        raise

    run_summary = output["summary"].get("run", {})
    run_messages = sum(int(v["messages"]) for v in run_summary.values())
    run_collectives = sum(int(v["collectives"]) for v in run_summary.values())
    run_bytes = sum(int(v["bytes"]) for v in run_summary.values())
    print(
        f"WSE plan merged: ranks={args.expected_ranks} "
        f"run_messages={run_messages} run_collectives={run_collectives} "
        f"run_bytes={run_bytes} output={output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
