#!/usr/bin/env python3
"""Estimate per-rank LAMMPS communication (T1) and computation (C1)."""

import argparse
import json
import math
import shlex
from pathlib import Path

QQRD2E_METAL = 14.3996454784255
SPHERE = 4.0 * math.pi / 3.0
PPPM_ACONS_5 = (
    1.0 / 23232.0,
    7601.0 / 13628160.0,
    143.0 / 69120.0,
    517231.0 / 106536960.0,
    106640677.0 / 11737571328.0,
)
EAM_PAIR_OPS_RANGE = (62.0, 66.0)


def commands(path):
    result = []
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        text = raw.split("#", 1)[0].strip()
        if text:
            result.append((lineno, shlex.split(text)))
    return result


def one(cmds, name, required=True):
    matches = [(line, args) for line, args in cmds if args[0] == name]
    if not matches:
        if required:
            raise ValueError(f"missing '{name}' command")
        return None
    return matches[-1]


def read_data(path, atom_style):
    lines = path.read_text().splitlines()
    count = None
    bounds = {}
    atoms = []
    for raw in lines:
        fields = raw.split()
        if len(fields) >= 2 and fields[1] == "atoms":
            count = int(fields[0])
        if len(fields) >= 4 and fields[-2:] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
            bounds[fields[-2][0]] = (float(fields[0]), float(fields[1]))

    start = next((i + 1 for i, text in enumerate(lines) if text.strip().startswith("Atoms")), None)
    if start is None:
        raise ValueError(f"{path}: no Atoms section")
    for raw in lines[start:]:
        text = raw.split("#", 1)[0].strip()
        if not text:
            continue
        fields = text.split()
        if not fields[0].lstrip("+-").isdigit():
            break
        if atom_style == "charge":
            atom_id, atom_type = int(fields[0]), int(fields[1])
            charge, xyz = float(fields[2]), tuple(map(float, fields[3:6]))
        elif atom_style == "atomic":
            atom_id, atom_type = int(fields[0]), int(fields[1])
            charge, xyz = 0.0, tuple(map(float, fields[2:5]))
        else:
            raise ValueError(f"unsupported atom_style: {atom_style}")
        atoms.append([atom_id, atom_type, charge, *xyz])
    if count != len(atoms) or set(bounds) != {"x", "y", "z"}:
        raise ValueError(f"{path}: incomplete atom data or orthogonal box")
    return atoms, [bounds[d][0] for d in "xyz"], [bounds[d][1] for d in "xyz"]


def make_fcc(cmds):
    _, lattice = one(cmds, "lattice")
    _, region = one(cmds, "region")
    if lattice[1] != "fcc" or region[1:3] != ["box", "block"]:
        raise ValueError("only 'lattice fcc' with 'region box block' is supported")
    a = float(lattice[2])
    limits = list(map(float, region[3:9]))
    lo = [limits[0] * a, limits[2] * a, limits[4] * a]
    hi = [limits[1] * a, limits[3] * a, limits[5] * a]
    cells = [round(limits[1] - limits[0]), round(limits[3] - limits[2]), round(limits[5] - limits[4])]
    if any(n <= 0 for n in cells):
        raise ValueError("FCC region must contain positive integer cell counts")
    basis = ((0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0))
    atoms = []
    for iz in range(cells[2]):
        for iy in range(cells[1]):
            for ix in range(cells[0]):
                for bx, by, bz in basis:
                    atoms.append([len(atoms) + 1, 1, 0.0,
                                  lo[0] + (ix + bx) * a,
                                  lo[1] + (iy + by) * a,
                                  lo[2] + (iz + bz) * a])
    return atoms, lo, hi


def replicate(atoms, lo, hi, factors):
    lengths = [hi[d] - lo[d] for d in range(3)]
    result = []
    for iz in range(factors[2]):
        for iy in range(factors[1]):
            for ix in range(factors[0]):
                shift = [ix * lengths[0], iy * lengths[1], iz * lengths[2]]
                for atom in atoms:
                    result.append([len(result) + 1, atom[1], atom[2],
                                   atom[3] + shift[0], atom[4] + shift[1], atom[5] + shift[2]])
    return result, lo, [lo[d] + lengths[d] * factors[d] for d in range(3)]


def apply_sets(cmds, atoms):
    groups = {}
    for _, args in cmds:
        if args[0] == "group" and len(args) >= 4 and args[2] == "id":
            ids = set()
            for spec in args[3:]:
                parts = list(map(int, spec.split(":")))
                if len(parts) == 1:
                    ids.add(parts[0])
                else:
                    start, stop = parts[:2]
                    stride = parts[2] if len(parts) == 3 else 1
                    ids.update(range(start, stop + 1, stride))
            groups[args[1]] = ids
        elif args[:2] == ["set", "group"] and len(args) == 5 and args[3] == "type":
            selected = groups.get(args[2], set())
            for atom in atoms:
                if atom[0] in selected:
                    atom[1] = int(args[4])
        elif args[:2] == ["set", "type"] and len(args) == 5 and args[3] == "charge":
            atom_type, charge = int(args[2]), float(args[4])
            for atom in atoms:
                if atom[1] == atom_type:
                    atom[2] = charge


def parse_system(input_path):
    cmds = commands(input_path)
    _, units = one(cmds, "units")
    if units[1] != "metal":
        raise ValueError("only 'units metal' is supported")
    _, boundary = one(cmds, "boundary")
    if boundary[1:4] != ["p", "p", "p"]:
        raise ValueError("only periodic 'boundary p p p' is supported")
    _, atom_style_cmd = one(cmds, "atom_style")
    atom_style = atom_style_cmd[1]
    data_cmd = one(cmds, "read_data", False)
    if data_cmd:
        data_path = (input_path.parent / data_cmd[1][1]).resolve()
        atoms, lo, hi = read_data(data_path, atom_style)
    else:
        atoms, lo, hi = make_fcc(cmds)
        data_path = None

    rep = one(cmds, "replicate", False)
    if rep:
        factors = tuple(map(int, rep[1][1:4]))
        atoms, lo, hi = replicate(atoms, lo, hi, factors)
    apply_sets(cmds, atoms)

    _, proc = one(cmds, "processors")
    procgrid = tuple(map(int, proc[1:4]))
    if math.prod(procgrid) <= 0:
        raise ValueError("invalid processors grid")
    _, pair = one(cmds, "pair_style")
    pair_style = pair[1]
    if pair_style == "eam":
        _, coeff = one(cmds, "pair_coeff")
        potential = (input_path.parent / coeff[-1]).resolve()
        header = potential.read_text().splitlines()[2].split()
        force_cutoff = float(header[4])
    elif pair_style.startswith("lj/cut/coul/"):
        force_cutoff = max(map(float, pair[2:]))
    else:
        raise ValueError(f"unsupported pair_style: {pair_style}")

    neighbor = one(cmds, "neighbor", False)
    skin = float(neighbor[1][1]) if neighbor else 2.0
    kspace = one(cmds, "kspace_style", False)
    if kspace and kspace[1][1] != "pppm":
        raise ValueError("only kspace_style pppm is supported")
    return {
        "commands": cmds,
        "atoms": atoms,
        "lo": lo,
        "hi": hi,
        "data_file": str(data_path) if data_path else None,
        "atom_style": atom_style,
        "procgrid": procgrid,
        "pair_style": pair_style,
        "force_cutoff": force_cutoff,
        "skin": skin,
        "kspace_accuracy": float(kspace[1][2]) if kspace else None,
    }


def rank_counts(system):
    lo, hi = system["lo"], system["hi"]
    lengths = [hi[d] - lo[d] for d in range(3)]
    px, py, pz = system["procgrid"]
    counts = [0] * (px * py * pz)
    for atom in system["atoms"]:
        loc = []
        for d, p in enumerate((px, py, pz)):
            cell = int((atom[d + 3] - lo[d]) / lengths[d] * p)
            loc.append(min(p - 1, max(0, cell)))
        rank = loc[0] * py * pz + loc[1] * pz + loc[2]
        counts[rank] += 1
    return counts


def expanded_copies(lengths, procgrid, width):
    local = [lengths[d] / procgrid[d] for d in range(3)]
    copies = 0.0
    details = []
    for d in range(3):
        if procgrid[d] == 1:
            continue
        rounds = int(width / local[d]) + 1
        area = math.prod(
            min(lengths[j], local[j] + 2.0 * width) if j < d else local[j]
            for j in range(3) if j != d
        )
        for round_id in range(rounds):
            slab = max(0.0, min(local[d], width - round_id * local[d]))
            if slab:
                volume = 2.0 * slab * area
                copies += volume
                details.append({"dimension": d, "round": round_id, "volume": volume})
    return copies, details


def commbrick_copies(system, width):
    """Replay CommBrick borders() selection and return network atom-copy counts."""
    lo, hi, procgrid = system["lo"], system["hi"], system["procgrid"]
    lengths = [hi[d] - lo[d] for d in range(3)]
    local = [lengths[d] / procgrid[d] for d in range(3)]
    ranks = math.prod(procgrid)
    owned = [[] for _ in range(ranks)]
    for atom in system["atoms"]:
        loc = [min(procgrid[d] - 1, max(0, int(
            (atom[d + 3] - lo[d]) / lengths[d] * procgrid[d]
        ))) for d in range(3)]
        rank = loc[0] * procgrid[1] * procgrid[2] + loc[1] * procgrid[2] + loc[2]
        owned[rank].append(tuple(atom[3:6]))
    present = [list(atoms) for atoms in owned]
    copies = [0] * ranks

    def location(rank):
        return [rank // (procgrid[1] * procgrid[2]),
                (rank // procgrid[2]) % procgrid[1],
                rank % procgrid[2]]

    for dim in range(3):
        maxneed = int(width * procgrid[dim] / lengths[dim]) + 1
        previous_end = [0] * ranks
        candidate_end = [0] * ranks
        for ineed in range(2 * maxneed):
            if ineed % 2 == 0:
                previous_end = candidate_end
                candidate_end = [len(atoms) for atoms in present]
            sends = [[] for _ in range(ranks)]
            destinations = [0] * ranks
            for rank in range(ranks):
                loc = location(rank)
                sublo = lo[dim] + loc[dim] * local[dim]
                subhi = sublo + local[dim]
                if ineed % 2 == 0:
                    slab_lo = -math.inf if ineed < 2 else 0.5 * (sublo + subhi)
                    slab_hi = sublo + width
                    dstloc = (loc[dim] - 1) % procgrid[dim]
                    pbc_shift = lengths[dim] if loc[dim] == 0 else 0.0
                else:
                    slab_lo = subhi - width
                    slab_hi = math.inf if ineed < 2 else 0.5 * (sublo + subhi)
                    dstloc = (loc[dim] + 1) % procgrid[dim]
                    pbc_shift = -lengths[dim] if loc[dim] == procgrid[dim] - 1 else 0.0
                dstcoords = loc[:]
                dstcoords[dim] = dstloc
                destinations[rank] = (dstcoords[0] * procgrid[1] * procgrid[2] +
                                      dstcoords[1] * procgrid[2] + dstcoords[2])
                for xyz in present[rank][previous_end[rank]:candidate_end[rank]]:
                    if slab_lo <= xyz[dim] <= slab_hi:
                        shifted = list(xyz)
                        shifted[dim] += pbc_shift
                        sends[rank].append(tuple(shifted))
            for rank, payload in enumerate(sends):
                destination = destinations[rank]
                present[destination].extend(payload)
                if destination != rank:
                    copies[rank] += len(payload)
    return copies


def pppm_grid(system, order=5):
    atoms = system["atoms"]
    natoms = len(atoms)
    lengths = [system["hi"][d] - system["lo"][d] for d in range(3)]
    accuracy = system["kspace_accuracy"] * QQRD2E_METAL
    qsqsum = sum(atom[2] * atom[2] for atom in atoms)
    q2 = qsqsum * QQRD2E_METAL
    cutoff = system["force_cutoff"]
    volume = math.prod(lengths)
    small = accuracy * math.sqrt(natoms * cutoff * volume) / (2.0 * q2)
    gewald = ((1.35 - 0.15 * math.log(accuracy)) if small >= 1.0
              else math.sqrt(-math.log(small))) / cutoff

    def error(h, period):
        hg = h * gewald
        poly = sum(PPPM_ACONS_5[m] * hg ** (2 * m) for m in range(order))
        return q2 * hg ** order * math.sqrt(
            gewald * period * math.sqrt(2.0 * math.pi) * poly / natoms
        ) / (period * period)

    grid = []
    for period in lengths:
        n = int(period * gewald) + 1
        h = 1.0 / gewald
        err = error(h, period)
        while err > accuracy:
            err = error(h, period)
            n += 1
            h = period / n
        while True:
            rest = n
            for factor in (2, 3, 5):
                while rest % factor == 0:
                    rest //= factor
            if rest == 1:
                break
            n += 1
        grid.append(n)
    return tuple(grid), gewald


def grid_halo_copies(grid, procgrid, order=5):
    local = [grid[d] / procgrid[d] for d in range(3)]
    halo = order / 2.0
    copies = 0.0
    for d in range(3):
        if procgrid[d] == 1:
            continue
        area = math.prod(
            min(grid[j], local[j] + 2.0 * halo) if j < d else local[j]
            for j in range(3) if j != d
        )
        copies += 2.0 * halo * area
    return copies


def estimate(system):
    atoms = system["atoms"]
    natoms = len(atoms)
    procgrid = system["procgrid"]
    ranks = math.prod(procgrid)
    lengths = [system["hi"][d] - system["lo"][d] for d in range(3)]
    volume = math.prod(lengths)
    density = natoms / volume
    local_volume = volume / ranks
    counts = rank_counts(system)
    ghost_width = system["force_cutoff"] + system["skin"]
    _, swaps = expanded_copies(lengths, procgrid, ghost_width)
    rank_copies = commbrick_copies(system, ghost_width)
    expanded_volume = math.prod(lengths[d] / procgrid[d] + 2.0 * ghost_width for d in range(3))
    average_ghosts = max(0.0, density * expanded_volume - natoms / ranks)

    kspace = None
    if system["kspace_accuracy"] is not None:
        grid, gewald = pppm_grid(system)
        mesh = math.prod(grid)
        grid_copies = grid_halo_copies(grid, procgrid)
        nonlocal_fraction = 1.0 - 1.0 / ranks
        kspace = {
            "grid": grid,
            "order": 5,
            "gewald": gewald,
            "mesh_points": mesh,
            "grid_reverse_bytes_per_rank": grid_copies * 8.0,
            "grid_forward_bytes_per_rank": grid_copies * 3.0 * 8.0,
            "fft_remap_bytes_per_rank": mesh / ranks * 160.0 * nonlocal_fraction,
        }

    rank_results = []
    for rank, nlocal in enumerate(counts):
        scale = nlocal / (natoms / ranks)
        copies = rank_copies[rank]
        nghost = average_ghosts * scale
        nlist = 0.5 * nlocal * density * SPHERE * ghost_width ** 3
        npair = 0.5 * nlocal * density * SPHERE * system["force_cutoff"] ** 3

        forward = copies * 3.0 * 8.0
        reverse = forward
        borders = copies * (7 if system["atom_style"] == "charge" else 6) * 8.0
        pair_comm = copies * 2.0 * 8.0 if system["pair_style"] == "eam" else 0.0
        comm = {
            "forward_bytes": forward,
            "reverse_bytes": reverse,
            "pair_forward_reverse_bytes": pair_comm,
            "borders_bytes": borders,
            "grid_reverse_bytes": 0.0,
            "grid_forward_bytes": 0.0,
            "fft_remap_bytes": 0.0,
        }

        integrate_ops = 18.0 * nlocal
        if system["pair_style"] == "eam":
            pair_min = (
                8.0 * nlist
                + EAM_PAIR_OPS_RANGE[0] * npair
                + 20.0 * (nlocal + nghost)
            )
            pair_max = (
                8.0 * nlist
                + EAM_PAIR_OPS_RANGE[1] * npair
                + 20.0 * (nlocal + nghost)
            )
        else:
            pair_min = 8.0 * nlist + 45.0 * npair
            pair_max = 8.0 * nlist + 50.0 * npair
        pppm_ops = 0.0
        if kspace:
            for key in ("grid_reverse_bytes_per_rank", "grid_forward_bytes_per_rank",
                        "fft_remap_bytes_per_rank"):
                comm[key.removesuffix("_per_rank")] = kspace[key] * scale
            mesh_local = kspace["mesh_points"] / ranks * scale
            particle_grid = nlocal * kspace["order"] ** 3
            pppm_ops = (40.0 * particle_grid +
                        20.0 * mesh_local * math.log2(kspace["mesh_points"]) +
                        20.0 * mesh_local)

        steady_t1 = (comm["forward_bytes"] + comm["reverse_bytes"] +
                     comm["pair_forward_reverse_bytes"] + comm["grid_reverse_bytes"] +
                     comm["grid_forward_bytes"] + comm["fft_remap_bytes"])
        rebuild_t1 = (comm["borders_bytes"] + comm["reverse_bytes"] +
                      comm["pair_forward_reverse_bytes"] + comm["grid_reverse_bytes"] +
                      comm["grid_forward_bytes"] + comm["fft_remap_bytes"])
        compute_min = integrate_ops + pair_min + pppm_ops
        compute_max = integrate_ops + pair_max + pppm_ops
        compute_midpoint = 0.5 * (compute_min + compute_max)
        rank_results.append({
            "rank": rank,
            "nlocal": nlocal,
            "nghost_estimate": nghost,
            "neighbor_entries_estimate": nlist,
            "force_pairs_estimate": npair,
            "communication": comm,
            "T1_send_bytes": steady_t1,
            "T1_steady_send_bytes": steady_t1,
            "T1_rebuild_send_bytes": rebuild_t1,
            "computation": {
                "integrate_ops": integrate_ops,
                "pair_ops_min": pair_min,
                "pair_ops_max": pair_max,
                "pppm_ops": pppm_ops,
            },
            "C1_ops": compute_midpoint,
            "C1_steady_ops_min": compute_min,
            "C1_steady_ops_max": compute_max,
            "C1_steady_ops_midpoint": compute_midpoint,
        })

    def stats(field):
        values = [rank[field] for rank in rank_results]
        return {"sum": sum(values), "average": sum(values) / ranks,
                "minimum": min(values), "maximum": max(values)}

    return {
        "schema_version": 1,
        "model": "static_uniform_density_v1",
        "input": {
            "atoms": natoms,
            "box": lengths,
            "density": density,
            "procgrid": procgrid,
            "num_ranks": ranks,
            "atom_style": system["atom_style"],
            "pair_style": system["pair_style"],
            "force_cutoff": system["force_cutoff"],
            "neighbor_skin": system["skin"],
            "ghost_width": ghost_width,
            "kspace_accuracy": system["kspace_accuracy"],
        },
        "assumptions": {
            "scope": "one steady Verlet timestep",
            "communication_metric": "application send payload; collectives excluded",
            "C1_metric": (
                "modeled double-precision arithmetic-equivalent operations; "
                "validate against lane-weighted FP_ARITH_INST_RETIRED"
            ),
            "newton_pair": True,
            "uniform_density_for_neighbors_and_ghosts": True,
            "pair_candidate_ops": 8,
            "lj_coul_kernel_ops_range": [45, 50],
            "eam_pair_ops_range": list(EAM_PAIR_OPS_RANGE),
            "pppm_fft_remap": "160 bytes/grid-point × (1-1/ranks), topology approximation",
            "commbrick_copies": "coordinate-level replay of CommBrick borders swaps",
        },
        "swap_geometry": swaps,
        "pppm": kspace,
        "summary": {
            "T1_send_bytes": stats("T1_send_bytes"),
            "T1_steady_send_bytes": stats("T1_steady_send_bytes"),
            "T1_rebuild_send_bytes": stats("T1_rebuild_send_bytes"),
            "C1_ops": stats("C1_ops"),
            "C1_steady_ops_min": stats("C1_steady_ops_min"),
            "C1_steady_ops_max": stats("C1_steady_ops_max"),
            "C1_steady_ops_midpoint": stats("C1_steady_ops_midpoint"),
        },
        "ranks": rank_results,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Estimate per-rank T1 send bytes and C1 scalar-equivalent operations")
    parser.add_argument("input", type=Path, help="LAMMPS input file")
    parser.add_argument("-o", "--output", type=Path, help="JSON output path")
    args = parser.parse_args()
    input_path = args.input.resolve()
    result = estimate(parse_system(input_path))
    result["input_file"] = str(input_path)
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")
    summary = result["summary"]
    print(
        f"T1 steady avg={summary['T1_steady_send_bytes']['average']:.0f} B/rank/step; "
        f"C1 avg={summary['C1_steady_ops_midpoint']['average']:.0f} ops/rank/step",
        file=__import__("sys").stderr,
    )


if __name__ == "__main__":
    main()
