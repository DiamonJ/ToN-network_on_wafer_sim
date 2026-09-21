#!/usr/bin/env python3
"""Derive uncalibrated per-rank LAMMPS T1/C1 theoretical bounds."""

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
# Source-counted scalar arithmetic conventions.  Each add/subtract, multiply,
# divide/reciprocal, sqrt, and exp is one algorithmic operation.  The bounds
# reflect branches visible in the kernels (table/direct Coulomb, energy/virial,
# and Newton updates); they are not fitted to a hardware profile.
PAIR_DISTANCE_OPS = {
    "coordinate_subtractions": 3,
    "squared_distance_multiplications": 3,
    "squared_distance_additions": 2,
}
EAM_DENSITY_OPS = {
    "radial_coordinate": 3,
    "two_cubic_interpolations_and_accumulations": 14,
}
EAM_FORCE_OPS = {
    "radial_coordinate": 3,
    "three_quadratic_interpolations": 12,
    "one_cubic_interpolation": 6,
    "reciprocal_and_pair_terms": 11,
    "force_accumulation": 12,
}
EAM_EMBED_OPS_RANGE = (6, 14)
EAM_ENERGY_VIRIAL_EXTRA_OPS = 12
LJ_COUL_FORCE_OPS_RANGE = (27, 44)
LJ_COUL_ENERGY_VIRIAL_EXTRA_OPS = 20
PPPM_PARTICLE_GRID_OPS_RANGE = (8, 20)
PPPM_FFT_BUTTERFLY_OPS_RANGE = (10, 20)
PPPM_MESH_FIELD_OPS_RANGE = (12, 20)
FFT_SCALAR_BYTES = 8
FFT_COMPLEX_SCALARS = 2
PPPM_FFT_TRANSFORMS_IK = 4  # one forward density FFT plus three inverse gradients
FFT_REMAPS_PER_TRANSFORM_RANGE = (2, 4)  # mid1/mid2; optional pre/post


def operation_sum(parts):
    return float(sum(parts.values()))


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
        kspace = {
            "grid": grid,
            "order": 5,
            "gewald": gewald,
            "mesh_points": mesh,
            "grid_reverse_bytes_per_rank": grid_copies * 8.0,
            "grid_forward_bytes_per_rank": grid_copies * 3.0 * 8.0,
            # Each FFT has two mandatory mid-remaps and up to two optional
            # pre/post remaps (fft3d.cpp).  Each complex mesh value carries two
            # FFT_SCALARs.  brick2fft adds one real-scalar remap.  The lower
            # remote-traffic bound is zero because a remap may be local for a
            # particular decomposition; the upper bound sends every element
            # for every possible remap to another rank.
            "fft_remap_bytes_lower_per_rank": 0.0,
            "fft_remap_bytes_upper_per_rank": (
                mesh / ranks
                * (
                    FFT_SCALAR_BYTES
                    + PPPM_FFT_TRANSFORMS_IK
                    * FFT_REMAPS_PER_TRANSFORM_RANGE[1]
                    * FFT_COMPLEX_SCALARS
                    * FFT_SCALAR_BYTES
                )
            ),
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
            "fft_remap_bytes_lower": 0.0,
            "fft_remap_bytes_upper": 0.0,
        }

        integrate_ops = 18.0 * nlocal
        distance_ops = operation_sum(PAIR_DISTANCE_OPS)
        if system["pair_style"] == "eam":
            # PairEAM traverses the neighbor list once for density and once for
            # forces, so the distance expression is evaluated twice.
            eam_pair_lower = operation_sum(EAM_DENSITY_OPS) + operation_sum(EAM_FORCE_OPS)
            eam_pair_upper = eam_pair_lower + EAM_ENERGY_VIRIAL_EXTRA_OPS
            pair_min = (
                2.0 * distance_ops * nlist
                + eam_pair_lower * npair
                + EAM_EMBED_OPS_RANGE[0] * nlocal
            )
            pair_max = (
                2.0 * distance_ops * nlist
                + eam_pair_upper * npair
                + EAM_EMBED_OPS_RANGE[1] * nlocal
            )
        else:
            pair_min = distance_ops * nlist + LJ_COUL_FORCE_OPS_RANGE[0] * npair
            pair_max = distance_ops * nlist + (
                LJ_COUL_FORCE_OPS_RANGE[1] + LJ_COUL_ENERGY_VIRIAL_EXTRA_OPS
            ) * npair
        pppm_min = 0.0
        pppm_max = 0.0
        if kspace:
            for key in ("grid_reverse_bytes_per_rank", "grid_forward_bytes_per_rank"):
                comm[key.removesuffix("_per_rank")] = kspace[key] * scale
            comm["fft_remap_bytes_lower"] = (
                kspace["fft_remap_bytes_lower_per_rank"] * scale
            )
            comm["fft_remap_bytes_upper"] = (
                kspace["fft_remap_bytes_upper_per_rank"] * scale
            )
            # Compatibility/budget field: use the conservative theoretical upper bound.
            comm["fft_remap_bytes"] = comm["fft_remap_bytes_upper"]
            mesh_local = kspace["mesh_points"] / ranks * scale
            particle_grid = nlocal * kspace["order"] ** 3
            fft_stages = math.log2(kspace["mesh_points"])
            pppm_min = (
                PPPM_PARTICLE_GRID_OPS_RANGE[0] * particle_grid
                + PPPM_FFT_BUTTERFLY_OPS_RANGE[0] * mesh_local * fft_stages
                + PPPM_MESH_FIELD_OPS_RANGE[0] * mesh_local
            )
            pppm_max = (
                PPPM_PARTICLE_GRID_OPS_RANGE[1] * particle_grid
                + PPPM_FFT_BUTTERFLY_OPS_RANGE[1] * mesh_local * fft_stages
                + PPPM_MESH_FIELD_OPS_RANGE[1] * mesh_local
            )

        steady_without_fft = (
            comm["forward_bytes"] + comm["reverse_bytes"]
            + comm["pair_forward_reverse_bytes"] + comm["grid_reverse_bytes"]
            + comm["grid_forward_bytes"]
        )
        rebuild_without_fft = (
            comm["borders_bytes"] + comm["reverse_bytes"]
            + comm["pair_forward_reverse_bytes"] + comm["grid_reverse_bytes"]
            + comm["grid_forward_bytes"]
        )
        steady_t1_lower = steady_without_fft + comm["fft_remap_bytes_lower"]
        steady_t1_upper = steady_without_fft + comm["fft_remap_bytes_upper"]
        rebuild_t1_lower = rebuild_without_fft + comm["fft_remap_bytes_lower"]
        rebuild_t1_upper = rebuild_without_fft + comm["fft_remap_bytes_upper"]
        # Compatibility/budget fields use the conservative upper endpoint.
        steady_t1 = steady_t1_upper
        rebuild_t1 = rebuild_t1_upper
        compute_min = integrate_ops + pair_min + pppm_min
        compute_max = integrate_ops + pair_max + pppm_max
        t1_lower = min(steady_t1_lower, rebuild_t1_lower)
        t1_upper = max(steady_t1_upper, rebuild_t1_upper)
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
            "T1_send_bytes_lower": t1_lower,
            "T1_send_bytes_upper": t1_upper,
            "T1_steady_send_bytes_lower": steady_t1_lower,
            "T1_steady_send_bytes_upper": steady_t1_upper,
            "T1_rebuild_send_bytes_lower": rebuild_t1_lower,
            "T1_rebuild_send_bytes_upper": rebuild_t1_upper,
            "computation": {
                "integrate_ops": integrate_ops,
                "pair_ops_min": pair_min,
                "pair_ops_max": pair_max,
                "pppm_ops_min": pppm_min,
                "pppm_ops_max": pppm_max,
            },
            "C1_steady_ops_min": compute_min,
            "C1_steady_ops_max": compute_max,
            "C1_steady_ops_lower": compute_min,
            "C1_steady_ops_upper": compute_max,
            # A simulator needs one scheduling budget.  Use the derived upper
            # bound, never a fitted midpoint or profile-derived coefficient.
            "C1_steady_ops_budget": compute_max,
        })

    def stats(field):
        values = [rank[field] for rank in rank_results]
        return {"sum": sum(values), "average": sum(values) / ranks,
                "minimum": min(values), "maximum": max(values)}

    return {
        "schema_version": 2,
        "model": "static_uncalibrated_bounds_v2",
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
            "communication_metric": (
                "application send payload bytes; collectives and atom migration "
                "exchange are excluded"
            ),
            "C1_metric": (
                "algorithmic scalar-equivalent operation bounds; not converted "
                "to retired instructions with a fitted coefficient"
            ),
            "calibration": "none",
            "bounds_policy": (
                "T1 lower/upper span steady versus neighbor-rebuild traffic; "
                "C1 lower/upper come from explicit kernel operation-count ranges"
            ),
            "newton_pair": True,
            "uniform_density_for_neighbors_and_ghosts": True,
            "operation_counting_convention": (
                "one scalar-equivalent operation per source add/subtract, "
                "multiply, divide/reciprocal, sqrt, or exp"
            ),
            "pair_distance_operation_breakdown": PAIR_DISTANCE_OPS,
            "lj_coul_force_ops_range": list(LJ_COUL_FORCE_OPS_RANGE),
            "lj_coul_energy_virial_extra_ops_upper": LJ_COUL_ENERGY_VIRIAL_EXTRA_OPS,
            "eam_density_operation_breakdown": EAM_DENSITY_OPS,
            "eam_force_operation_breakdown": EAM_FORCE_OPS,
            "eam_embedding_ops_range_per_atom": list(EAM_EMBED_OPS_RANGE),
            "eam_energy_virial_extra_ops_upper": EAM_ENERGY_VIRIAL_EXTRA_OPS,
            "pppm_particle_grid_ops_range": list(PPPM_PARTICLE_GRID_OPS_RANGE),
            "pppm_fft_butterfly_ops_range": list(PPPM_FFT_BUTTERFLY_OPS_RANGE),
            "pppm_mesh_field_ops_range": list(PPPM_MESH_FIELD_OPS_RANGE),
            "pppm_fft_remap": {
                "fft_scalar_bytes": FFT_SCALAR_BYTES,
                "complex_scalars": FFT_COMPLEX_SCALARS,
                "transforms_ik": PPPM_FFT_TRANSFORMS_IK,
                "remaps_per_transform_range": list(FFT_REMAPS_PER_TRANSFORM_RANGE),
                "remote_fraction_range": [0.0, 1.0],
            },
            "commbrick_copies": "coordinate-level replay of CommBrick borders swaps",
        },
        "swap_geometry": swaps,
        "pppm": kspace,
        "summary": {
            "T1_send_bytes": stats("T1_send_bytes"),
            "T1_steady_send_bytes": stats("T1_steady_send_bytes"),
            "T1_rebuild_send_bytes": stats("T1_rebuild_send_bytes"),
            "T1_send_bytes_lower": stats("T1_send_bytes_lower"),
            "T1_send_bytes_upper": stats("T1_send_bytes_upper"),
            "T1_steady_send_bytes_lower": stats("T1_steady_send_bytes_lower"),
            "T1_steady_send_bytes_upper": stats("T1_steady_send_bytes_upper"),
            "T1_rebuild_send_bytes_lower": stats("T1_rebuild_send_bytes_lower"),
            "T1_rebuild_send_bytes_upper": stats("T1_rebuild_send_bytes_upper"),
            "C1_steady_ops_min": stats("C1_steady_ops_min"),
            "C1_steady_ops_max": stats("C1_steady_ops_max"),
            "C1_steady_ops_lower": stats("C1_steady_ops_lower"),
            "C1_steady_ops_upper": stats("C1_steady_ops_upper"),
            "C1_steady_ops_budget": stats("C1_steady_ops_budget"),
        },
        "ranks": rank_results,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Derive per-rank T1 byte and C1 operation bounds without calibration")
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
        f"T1 avg bound=[{summary['T1_send_bytes_lower']['average']:.0f}, "
        f"{summary['T1_send_bytes_upper']['average']:.0f}] B/rank/step; "
        f"C1 avg bound=[{summary['C1_steady_ops_lower']['average']:.0f}, "
        f"{summary['C1_steady_ops_upper']['average']:.0f}] ops/rank/step",
        file=__import__("sys").stderr,
    )


if __name__ == "__main__":
    main()
