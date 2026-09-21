#!/usr/bin/env python3
"""Phase-1 compiler from LAMMPS WSE plans to deterministic WSE wavefronts.

The compiler deliberately accepts source-level communication plans, not MPI
traces.  Its first version covers the paper's short-range CommBrick domain:
2-D worker placement, fold-PBC, H/V stages, (b+1) strip phases, multicast
forest footprints, reverse reductions, link reservations, and EST output.
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import re
from pathlib import Path
from typing import Any


REVERSE_PHASES = frozenset(("reverse", "pair_reverse", "grid_reverse"))
SHORT_PHASES = frozenset(("forward", "reverse", "pair_forward", "pair_reverse"))


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as fh:
        value = json.load(fh)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON value must be an object")
    return value


def read_cfg(path: Path) -> dict[str, Any]:
    """Read BookSim's simple ``key = value;`` configuration syntax."""
    cfg: dict[str, Any] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("//", 1)[0].strip()
        match = re.match(r"^(?:\s*\d+\|)?\s*([A-Za-z_]\w*)\s*=\s*(.*?)\s*;\s*$", line)
        if not match:
            continue
        key, text = match.groups()
        if re.fullmatch(r"[-+]?\d+", text):
            value: Any = int(text)
        else:
            try:
                value = float(text)
            except ValueError:
                value = text
        cfg[key] = value
    required = ("topology", "k", "n", "flit_size_bytes", "noc_frequency_ghz")
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"{path}: missing configuration keys: {', '.join(missing)}")
    if cfg["topology"] != "mesh" or int(cfg["n"]) != 2:
        raise ValueError("Phase-1 WSE compiler requires topology=mesh and n=2")
    return cfg


def fold_coord(coord: int, size: int) -> int:
    """Paper III-E interleaving: periodic neighbors become <=2 mesh hops."""
    return 2 * coord if coord < size // 2 else 2 * (size - 1 - coord) + 1


def route_xy(src: tuple[int, int], dst: tuple[int, int]) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    links = []
    x, y = src
    while x != dst[0]:
        nxt = (x + (1 if dst[0] > x else -1), y)
        links.append(((x, y), nxt))
        x, y = nxt
    while y != dst[1]:
        nxt = (x, y + (1 if dst[1] > y else -1))
        links.append(((x, y), nxt))
        x, y = nxt
    return links


def phase_mode(phase: str) -> str:
    return "reduction" if phase in REVERSE_PHASES else "multicast"


class TreeLinkTable:
    """Reserve variable-length slots on every directed tree edge and ejection."""

    def __init__(self, hop_stride: int):
        self.hop_stride = hop_stride
        self.intervals: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
        self.retries = 0

    @staticmethod
    def _conflict(intervals: list[tuple[int, int, str]], start: int, end: int) -> int | None:
        escape = None
        for old_start, old_end, _ in intervals:
            if start <= old_end and end >= old_start:
                escape = max(escape or 0, old_end + 1)
        return escape

    def first_free(self, slots: list[dict[str, Any]], earliest: int) -> tuple[int, int]:
        candidate = earliest
        while True:
            escaped = candidate
            for slot in slots:
                offset = int(slot["depth"]) * self.hop_stride
                start = candidate + offset
                end = start + int(slot["flits"]) - 1
                conflict_end = self._conflict(self.intervals[slot["link"]], start, end)
                if conflict_end is not None:
                    escaped = max(escaped, conflict_end - offset)
            if escaped == candidate:
                completion = candidate + max(
                    int(slot["depth"]) * self.hop_stride + int(slot["flits"])
                    for slot in slots
                )
                return candidate, completion
            candidate = escaped
            self.retries += 1
            if candidate - earliest > (1 << 32):
                raise RuntimeError("TreeLinkTable search did not converge")

    def reserve(self, slots: list[dict[str, Any]], start: int, owner: str) -> None:
        for slot in slots:
            begin = start + int(slot["depth"]) * self.hop_stride
            end = begin + int(slot["flits"]) - 1
            self.intervals[slot["link"]].append((begin, end, owner))

    def assert_conflict_free(self) -> None:
        for link, intervals in self.intervals.items():
            ordered = sorted(intervals)
            for left, right in zip(ordered, ordered[1:]):
                if right[0] <= left[1]:
                    raise AssertionError(f"link conflict on {link}: {left} vs {right}")


def normalize_inputs(
    plan: dict[str, Any],
    cost: dict[str, Any],
    cfg: dict[str, Any],
    fold_pbc: bool = True,
) -> tuple[dict[str, Any], list[str]]:
    warnings = []
    ranks = int(plan.get("num_ranks", 0))
    procgrid = tuple(int(v) for v in plan.get("procgrid", ()))
    if len(procgrid) != 3 or math.prod(procgrid) != ranks:
        raise ValueError(f"invalid procgrid {procgrid} for {ranks} ranks")
    if procgrid[2] != 1:
        raise ValueError(f"Phase-1 requires Z compression (Pz=1), got {procgrid}")
    nx, ny, _ = procgrid
    if nx != ny or nx * ny != ranks:
        raise ValueError(f"Phase-1 BookSim target requires a square worker grid, got {nx}x{ny}")
    if int(cfg["k"]) != nx:
        warnings.append(
            f"cfg k={cfg['k']} is a template value; effective k={nx} comes from wse_plan procgrid"
        )
    fold = fold_pbc and nx % 2 == 0 and ny % 2 == 0
    if fold_pbc and not fold:
        warnings.append("fold-PBC disabled because both worker-grid dimensions must be even")

    metadata = plan.get("rank_metadata", [])
    if len(metadata) != ranks:
        raise ValueError(f"rank_metadata has {len(metadata)} entries, expected {ranks}")
    logical: dict[int, tuple[int, int]] = {}
    placement = []
    occupied = set()
    for item in metadata:
        rank = int(item["rank"])
        loc = tuple(int(v) for v in item["myloc"])
        if len(loc) != 3 or loc[2] != 0:
            raise ValueError(f"rank {rank}: invalid compressed coordinate {loc}")
        physical = (
            fold_coord(loc[0], nx) if fold else loc[0],
            fold_coord(loc[1], ny) if fold else loc[1],
        )
        if physical in occupied:
            raise AssertionError(f"placement is not bijective at {physical}")
        occupied.add(physical)
        logical[rank] = (loc[0], loc[1])
        placement.append(
            {"rank": rank, "logical": [loc[0], loc[1]], "physical": list(physical)}
        )
    if set(logical) != set(range(ranks)):
        raise ValueError("rank_metadata does not cover ranks [0,num_ranks)")

    cost_ranks = cost.get("ranks", [])
    if len(cost_ranks) != ranks:
        raise ValueError(f"cost model has {len(cost_ranks)} ranks, expected {ranks}")
    config = {
        "mesh": [nx, ny],
        "configured_k": int(cfg["k"]),
        "fold_pbc": fold,
        "flit_size_bytes": int(cfg["flit_size_bytes"]),
        "noc_frequency_ghz": float(cfg["noc_frequency_ghz"]),
        "compute_capability_ops_s": float(cfg.get("ccdg_compute_capability", 0)),
        "hop_stride_cycles": int(cfg.get("routing_delay", 0))
        + int(cfg.get("vc_alloc_delay", 0))
        + int(cfg.get("sw_alloc_delay", 0))
        + 1,
        "logical": logical,
        "placement": placement,
    }
    if config["flit_size_bytes"] <= 0:
        raise ValueError("flit_size_bytes must be positive")
    return config, warnings


def make_flow(
    record: dict[str, Any], config: dict[str, Any], placement_by_rank: dict[int, tuple[int, int]]
) -> dict[str, Any]:
    src, dst = int(record["src"]), int(record["dst"])
    path = route_xy(placement_by_rank[src], placement_by_rank[dst])
    if not path:
        raise ValueError(f"self-message is outside Phase-1 WSE semantics: rank {src}")
    if len({edge[0] for edge in path}) != len(path):
        raise AssertionError(f"cyclic route generated for {src}->{dst}")
    return {
        "src_rank": src,
        "dst_rank": dst,
        "src": list(placement_by_rank[src]),
        "dst": list(placement_by_rank[dst]),
        "round": int(record.get("round", 0)),
        "bytes": int(record.get("bytes", 0)),
        "path": [[list(a), list(b)] for a, b in path],
    }


def build_slots(flows: list[dict[str, Any]], flit_bytes: int) -> list[dict[str, Any]]:
    """Conservative aggregation: absent payload lineage, coincident loads sum."""
    loads: collections.Counter[tuple[str, int]] = collections.Counter()
    for flow in flows:
        src = flow["src"]
        loads[(f"INJECT@{src[0]},{src[1]}", 0)] += flow["bytes"]
        for depth, (src, dst) in enumerate(flow["path"]):
            loads[(f"{src[0]},{src[1]}->{dst[0]},{dst[1]}", depth)] += flow["bytes"]
        dst = flow["dst"]
        loads[(f"EJECT@{dst[0]},{dst[1]}", len(flow["path"]))] += flow["bytes"]
    slots = []
    for (link, depth), size in sorted(loads.items()):
        slots.append(
            {
                "link": link,
                "depth": depth,
                "bytes": size,
                "flits": max(1, math.ceil(size / flit_bytes)),
            }
        )
    return slots


def reduction_nodes(flows: list[dict[str, Any]]) -> list[list[int]]:
    indegree: collections.Counter[tuple[int, int]] = collections.Counter()
    outdegree: collections.Counter[tuple[int, int]] = collections.Counter()
    for flow in flows:
        for src, dst in flow["path"]:
            indegree[tuple(dst)] += 1
            outdegree[tuple(src)] += 1
    return [list(node) for node in sorted(indegree) if indegree[node] >= 2 and outdegree[node] >= 1]


def compile_program(
    plan_path: Path,
    cost_path: Path,
    cfg_path: Path,
    scope: str = "run",
    wse_fast_profile: bool = False,
    fold_pbc: bool = True,
    compute_capability: float | None = None,
) -> tuple[dict[str, Any], list[tuple[str, int]], dict[str, Any]]:
    plan, cost, cfg = read_json(plan_path), read_json(cost_path), read_cfg(cfg_path)
    if compute_capability is not None:
        if compute_capability <= 0:
            raise ValueError("compute capability override must be positive")
        cfg["ccdg_compute_capability"] = compute_capability
    if wse_fast_profile:
        cfg["flit_size_bytes"] = 4
        cfg["routing_delay"] = 0
        cfg["vc_alloc_delay"] = 1
        cfg["sw_alloc_delay"] = 1
    config, warnings = normalize_inputs(plan, cost, cfg, fold_pbc)
    config["profile"] = "wse_fast" if wse_fast_profile else "cfg_native"
    placement = {p["rank"]: tuple(p["physical"]) for p in config["placement"]}
    logical = config.pop("logical")

    all_messages = [
        item
        for item in plan.get("records", [])
        if item.get("kind") == "message" and item.get("scope") == scope
    ]
    unsupported = sorted({item.get("phase", "unknown") for item in all_messages} - SHORT_PHASES)
    if unsupported:
        raise ValueError(
            "Phase-1 short-range compiler cannot lower phases: " + ", ".join(unsupported)
        )
    messages = [
        item for item in all_messages if int(item.get("bytes", 0)) > 0 and int(item["src"]) != int(item["dst"])
    ]
    if not messages:
        raise ValueError(f"no non-empty {scope!r} short-range messages in {plan_path}")

    phase_first_seq: dict[tuple[int, str], int] = {}
    for item in messages:
        key = (int(item.get("timestep", 0)), str(item["phase"]))
        phase_first_seq[key] = min(phase_first_seq.get(key, 1 << 60), int(item.get("seq", 0)))

    rate = config["compute_capability_ops_s"] / (config["noc_frequency_ghz"] * 1e9)
    if rate <= 0:
        rate = float(cfg.get("ccdg_compute_rate", 1.0))
        warnings.append("compute capability is zero; using legacy ccdg_compute_rate")
    compute_blocks = []
    for item in cost["ranks"]:
        ops = float(item["C1_steady_ops_midpoint"])
        compute_blocks.append(
            {"rank": int(item["rank"]), "ops": ops, "cycles": math.ceil(ops / rate)}
        )
    compute_barrier = max(item["cycles"] for item in compute_blocks)

    table = TreeLinkTable(config["hop_stride_cycles"])
    timesteps = []
    est: list[tuple[str, int]] = []
    stage_serial = 0
    comm_open = compute_barrier
    for timestep in sorted({int(item.get("timestep", 0)) for item in messages}):
        timestep_stages = []
        phases = sorted(
            {str(item["phase"]) for item in messages if int(item.get("timestep", 0)) == timestep},
            key=lambda phase: (phase_first_seq[(timestep, phase)], phase),
        )
        for phase in phases:
            phase_messages = [
                item
                for item in messages
                if int(item.get("timestep", 0)) == timestep and item["phase"] == phase
            ]
            dimensions = sorted({int(item["dimension"]) for item in phase_messages})
            for dimension in dimensions:
                axis = "H" if dimension == 0 else "V"
                stage_messages = [item for item in phase_messages if int(item["dimension"]) == dimension]
                b = max(int(item.get("round", 0)) for item in stage_messages) + 1
                width = b + 1
                groups: dict[tuple[int, int], list[dict[str, Any]]] = collections.defaultdict(list)
                for item in stage_messages:
                    src_logical = logical[int(item["src"])]
                    perpendicular = src_logical[1] if dimension == 0 else src_logical[0]
                    groups[(perpendicular % width, int(item.get("direction", 0)))].append(item)

                stage_id = f"t{timestep}.{phase}.{axis}.{stage_serial}"
                stage_serial += 1
                wavefronts = []
                stage_open = comm_open
                for tick in range(width):
                    tick_open = stage_open
                    tick_finish = tick_open
                    for direction in sorted(direction for q, direction in groups if q == tick):
                        records = groups[(tick, direction)]
                        flows = [make_flow(item, config, placement) for item in records]
                        slots = build_slots(flows, config["flit_size_bytes"])
                        wave_id = f"{stage_id}.q{tick}.d{direction:+d}"
                        start, completion = table.first_free(slots, tick_open)
                        table.reserve(slots, start, wave_id)
                        wave = {
                            "id": wave_id,
                            "phase_tick": tick,
                            "direction": direction,
                            "start_cycle": start,
                            "completion_cycle": completion,
                            "flows": flows,
                            "footprint_slots": slots,
                            "logical_bytes": sum(flow["bytes"] for flow in flows),
                        }
                        if phase_mode(phase) == "reduction":
                            wave["reduction_op"] = "sum"
                            wave["reduction_nodes"] = reduction_nodes(flows)
                        wavefronts.append(wave)
                        est.append((wave_id, start))
                        tick_finish = max(tick_finish, completion)
                    stage_open = tick_finish
                stage = {
                    "id": stage_id,
                    "phase": phase,
                    "axis": axis,
                    "dimension": dimension,
                    "mode": phase_mode(phase),
                    "b": b,
                    "phase_count": width,
                    "start_cycle": comm_open,
                    "completion_cycle": stage_open,
                    "wavefronts": wavefronts,
                }
                timestep_stages.append(stage)
                comm_open = stage_open
        timesteps.append({"timestep": timestep, "stages": timestep_stages})

    table.assert_conflict_free()
    input_bytes = sum(int(item["bytes"]) for item in messages)
    emitted_bytes = sum(
        wave["logical_bytes"]
        for timestep in timesteps
        for stage in timestep["stages"]
        for wave in stage["wavefronts"]
    )
    if input_bytes != emitted_bytes:
        raise AssertionError(f"message conservation failed: input={input_bytes}, emitted={emitted_bytes}")
    est.sort(key=lambda item: (item[1], item[0]))
    starts = [start for _, start in est]
    if any(right < left for left, right in zip(starts, starts[1:])):
        raise AssertionError("EST is not monotonic")

    program = {
        "schema_version": 1,
        "compiler": "wse_phase1_short_v1",
        "inputs": {
            "wse_plan": str(plan_path),
            "static_cost": str(cost_path),
            "booksim_cfg": str(cfg_path),
            "scope": scope,
            "hardware_profile": config["profile"],
            "fold_pbc": fold_pbc,
            "compute_capability_override": compute_capability,
        },
        "hardware": config,
        "placement": config.pop("placement"),
        "compute_model": {
            "ops_per_noc_cycle": rate,
            "serialization": "global compute barrier before compiled communication",
            "blocks": compute_blocks,
        },
        "timesteps": timesteps,
    }
    max_fold_hops = max(
        len(flow["path"])
        for timestep in timesteps
        for stage in timestep["stages"]
        for wave in stage["wavefronts"]
        for flow in wave["flows"]
    )
    report = {
        "schema_version": 1,
        "status": "PASS",
        "scope": scope,
        "warnings": warnings,
        "limitations": [
            "source plan has aggregate message sizes but no payload lineage; coincident link loads are summed conservatively",
            "compute is a static global barrier; source-level compute/communication overlap is not yet emitted",
            "FFT remap, Grid3d, and collectives are outside the Phase-1 short-range compiler",
        ],
        "counts": {
            "ranks": int(plan["num_ranks"]),
            "messages": len(messages),
            "stages": sum(len(item["stages"]) for item in timesteps),
            "wavefronts": len(est),
            "link_reservation_retries": table.retries,
        },
        "conservation": {
            "input_messages": len(messages),
            "emitted_flows": sum(
                len(wave["flows"])
                for item in timesteps
                for stage in item["stages"]
                for wave in stage["wavefronts"]
            ),
            "input_bytes": input_bytes,
            "emitted_logical_bytes": emitted_bytes,
        },
        "cycles": {
            "compute_barrier": compute_barrier,
            "communication_completion": comm_open,
            "compiled_total": comm_open,
        },
        "assertions": {
            "placement_bijective": True,
            "z_compressed": True,
            "individual_routes_acyclic": True,
            "link_schedule_conflict_free": True,
            "est_monotonic": True,
            "message_count_conserved": len(messages)
            == sum(
                len(wave["flows"])
                for item in timesteps
                for stage in item["stages"]
                for wave in stage["wavefronts"]
            ),
            "message_bytes_conserved": input_bytes == emitted_bytes,
            "fold_neighbor_hops_le_2": (not config["fold_pbc"]) or max_fold_hops <= 2,
        },
    }
    if not all(report["assertions"].values()):
        raise AssertionError(f"compiler assertion failed: {report['assertions']}")
    return program, est, report


def build_replay_ccdg(
    program: dict[str, Any]
) -> tuple[dict[str, Any], list[tuple[int, int]], dict[str, int]]:
    """Expand compiler footprints into one-hop manager-level branch packets."""
    nx, ny = program["hardware"]["mesh"]
    rank_count = nx * ny
    logical_to_physical = {
        int(item["rank"]): int(item["physical"][0]) + nx * int(item["physical"][1])
        for item in program["placement"]
    }
    compute_by_logical = {
        int(item["rank"]): item for item in program["compute_model"]["blocks"]
    }
    physical_to_logical = {physical: logical for logical, physical in logical_to_physical.items()}
    nodes_by_rank: dict[int, list[tuple[int, dict[str, Any]]]] = collections.defaultdict(list)
    replay_events: list[tuple[int, dict[str, Any], int]] = []
    resource_intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    replay_stats = {
        "packets": 0,
        "commands": 0,
        "branches": 0,
        "bytes": 0,
        "branch_link_bytes": 0,
        "wavefronts": 0,
        "completion_cycle": 0,
    }

    for physical_rank in range(rank_count):
        logical_rank = physical_to_logical[physical_rank]
        block = compute_by_logical[logical_rank]
        nodes_by_rank[physical_rank].append(
            (
                0,
                {
                    "rank": physical_rank,
                    "type": "COMPUTE",
                    "compute_ops": block["ops"],
                    "compute_cycles": block["cycles"],
                },
            )
        )
    stage_idx = 0
    for timestep in program["timesteps"]:
        for stage in timestep["stages"]:
            for wave in stage["wavefronts"]:
                wave_idx = replay_stats["wavefronts"]
                replay_stats["wavefronts"] += 1
                release = int(wave["start_cycle"])
                root = logical_to_physical[int(wave["flows"][0]["src_rank"])]
                command_bytes = 0
                replay_events.append(
                    (
                        release,
                        {
                            "rank": root,
                            "type": "SEND",
                            "comm_src": root,
                            "comm_dst": root,
                            "comm_tag": replay_stats["packets"],
                            "comm_bytes": command_bytes,
                            "wse_wavefront_idx": wave_idx,
                            "wse_stage_idx": stage_idx,
                            "wse_kind": "command",
                            "wse_mode": stage["mode"],
                        },
                        1,
                    )
                )
                replay_stats["packets"] += 1
                replay_stats["commands"] += 1

                for slot in wave["footprint_slots"]:
                    match = re.fullmatch(
                        r"(\d+),(\d+)->(\d+),(\d+)", str(slot["link"])
                    )
                    if match is None:
                        continue
                    sx, sy, dx, dy = (int(value) for value in match.groups())
                    src, dst = sy * nx + sx, dy * nx + dx
                    size = int(slot["bytes"])
                    branch_release = release + int(slot["depth"]) * int(
                        program["hardware"]["hop_stride_cycles"]
                    )
                    replay_events.append(
                        (
                            branch_release,
                            {
                                "rank": src,
                                "type": "SEND",
                                "comm_src": src,
                                "comm_dst": dst,
                                "comm_tag": replay_stats["packets"],
                                "comm_bytes": size,
                                "wse_wavefront_idx": wave_idx,
                                "wse_stage_idx": stage_idx,
                                "wse_kind": "branch",
                                "wse_mode": stage["mode"],
                            },
                            int(slot["flits"]),
                        )
                    )
                    replay_stats["packets"] += 1
                    replay_stats["branches"] += 1
                    replay_stats["bytes"] += size
                    replay_stats["branch_link_bytes"] += size
            stage_idx += 1

    hop_stride = int(program["hardware"]["hop_stride_cycles"])
    for nominal, node, flits in sorted(
        replay_events,
        key=lambda event: (
            event[0],
            int(event[1]["wse_wavefront_idx"]),
            event[1]["wse_kind"],
            int(event[1]["rank"]),
            int(event[1]["comm_dst"]),
        ),
    ):
        src, dst = int(node["rank"]), int(node["comm_dst"])
        resources = [(f"INJECT@{src}", 0), (f"EJECT@{dst}", hop_stride)]
        if src != dst:
            resources.append((f"LINK@{src}->{dst}", hop_stride))
        candidate = nominal
        while True:
            escaped = candidate
            for resource, offset in resources:
                begin, end = candidate + offset, candidate + offset + flits - 1
                for old_begin, old_end in resource_intervals[resource]:
                    if begin <= old_end and end >= old_begin:
                        escaped = max(escaped, old_end - offset + 1)
            if escaped == candidate:
                break
            candidate = escaped
        for resource, offset in resources:
            resource_intervals[resource].append(
                (candidate + offset, candidate + offset + flits - 1)
            )
        nodes_by_rank[src].append((candidate, node))
        replay_stats["completion_cycle"] = max(
            replay_stats["completion_cycle"], candidate + hop_stride + flits
        )

    nodes = []
    schedule = []
    node_id = 0
    for rank in range(rank_count):
        ordered = sorted(
            enumerate(nodes_by_rank[rank]), key=lambda item: (item[1][0], item[0])
        )
        previous = None
        for _, (release, node) in ordered:
            node["id"] = node_id
            node["predecessors"] = [] if previous is None else [previous]
            nodes.append(node)
            schedule.append((node_id, release))
            previous = node_id
            node_id += 1
    replay = {
        "schema_version": 1,
        "format": "wse_manager_footprint_replay_v1",
        "num_ranks": rank_count,
        "wse_wavefronts": replay_stats["wavefronts"],
        "wse_commands": replay_stats["commands"],
        "wse_branches": replay_stats["branches"],
        "nodes": nodes,
        "cross_rank_edges": [],
    }
    return replay, schedule, replay_stats


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wse_plan", type=Path)
    parser.add_argument("static_cost", type=Path)
    parser.add_argument("booksim_cfg", type=Path)
    parser.add_argument("-o", "--output-prefix", type=Path, required=True)
    parser.add_argument("--scope", default="run", choices=("run", "setup"))
    parser.add_argument(
        "--wse-fast-profile",
        action="store_true",
        help="override the base cfg with 4-byte flits and the shallowest valid IQ pipeline",
    )
    parser.add_argument(
        "--no-fold",
        action="store_true",
        help="place logical workers directly instead of applying fold-PBC",
    )
    parser.add_argument(
        "--compute-capability",
        type=float,
        help="override ccdg_compute_capability from the base cfg (ops/s)",
    )
    args = parser.parse_args()

    program, est, report = compile_program(
        args.wse_plan,
        args.static_cost,
        args.booksim_cfg,
        args.scope,
        args.wse_fast_profile,
        not args.no_fold,
        args.compute_capability,
    )
    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    replay, replay_schedule, replay_stats = build_replay_ccdg(program)
    program_path = prefix.with_suffix(".program.json")
    est_path = prefix.with_suffix(".est")
    report_path = prefix.with_suffix(".report.json")
    replay_path = prefix.with_suffix(".replay.ccdg")
    replay_est_path = prefix.with_suffix(".replay.est")
    report["replay"] = replay_stats
    semantic_total = int(report["cycles"]["compiled_total"])
    replay_total = max(semantic_total, int(replay_stats["completion_cycle"]))
    report["cycles"]["wse_semantic_total"] = semantic_total
    report["cycles"]["manager_replay_total"] = replay_total
    report["cycles"]["compiled_total"] = replay_total
    program_path.write_text(json.dumps(program, indent=2) + "\n", encoding="utf-8")
    est_path.write_text(
        "# wavefront_id earliest_start_cycle\n"
        + "".join(f"{wave_id} {start}\n" for wave_id, start in est),
        encoding="utf-8",
    )
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    replay_path.write_text(json.dumps(replay, indent=2) + "\n", encoding="utf-8")
    replay_est_path.write_text(
        "# node_id earliest_start_cycle\n"
        + "".join(f"{node_id} {start}\n" for node_id, start in replay_schedule),
        encoding="utf-8",
    )
    print(
        "WSE compile PASS: "
        f"messages={report['counts']['messages']} stages={report['counts']['stages']} "
        f"wavefronts={report['counts']['wavefronts']} cycles={report['cycles']['compiled_total']}"
    )
    print(
        f"program={program_path} est={est_path} report={report_path} "
        f"replay={replay_path} replay_est={replay_est_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
