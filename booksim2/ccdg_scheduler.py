#!/usr/bin/env python3
"""CCDG -> deterministic schedule (EST) orchestrator.

Compute a compile-time Earliest-Start-Time schedule for every node in a
CCDG, using the same cost model as CCDGTrafficManager:

  - COMPUTE dwell   = compute_ops / rate      (rate = cap/(noc_ghz*1e9))
                      fallback compute_cycles / freq_ratio
  - collective dwell = expanded critical path (recursive doubling / binomial
                      tree rounds x (per-packet serialization + hop delay))
  - cross-edge flight = hops(src,dst) x hop_lat + flits(bytes) x ser_lat
                        (flit_size_bytes=1 by default, matching run_ccdg_mesh.sh)

Graph: each node's chain predecessors (same-rank `predecessors` array) plus
the cross_rank_edges (SEND -> WAITALL). EST is the longest-path value atas any
predecessor chain, i.e. the earliest cycle the node may legally start.

Outputs:
  <prefix>_sched.est    plain table "node_id est_cycles" (fed to BookSim)
  <prefix>_sorted.ccdg  full CCDG with nodes reordered by EST + sched_est field
  stdout report         critical-path makespan bound, predicted sched slots,
                        edge-delay distribution

Semantics note: an EST is a RELEASE time, not a deadline. If the simulation
network needs longer than the modeled flight delay, the PE falls back to
blocked (safe); if the schedule is conservative, the PE sits in sched_wait.
The orchestrator's accuracy only shifts the ledger between the two, never
breaks correctness.
"""
import argparse
import collections
import json
import math
import sys

COLLECTIVE_ROUNDS = {
    "ALLREDUCE": "rd", "ALLGATHER": "rd",
    "BCAST": "tree", "SCATTER": "tree",
    "REDUCE": "tree", "GATHER": "tree",
    "ALLTOALL": "alltoall",
}


def hop(a, b, k):
    """Manhattan distance on a k x k mesh, rank->node (r%k, r//k)."""
    return abs((a % k) - (b % k)) + abs((a // k) - (b // k))


def avg_hop(k):
    if k <= 1:
        return 0.0
    return 2.0 * (k * k - 1) / (3.0 * k)


def flits_of(bytes_, flit_size):
    return max(1, (bytes_ + flit_size - 1) // flit_size)


def collective_dur(node, N, k, ser_lat, hop_lat):
    """Expanded critical-path estimate for a collective node.

    comm_bytes carries the N-fold total; each expanded packet carries
    bytes/N (same division as CCDGTrafficManager::_advancePE).
    """
    t = node["type"]
    base = node.get("comm_bytes", 0) or 64
    per = max(1, base // N)
    ser = per * ser_lat

    if COLLECTIVE_ROUNDS.get(t) == "rd":
        # recursive doubling: round i pairs rank r with r^mask, additive
        # (reduction semantics -> rounds serialize). Take the max hop over
        # all ranks for the critical path.
        p = 1
        while p * 2 <= N:
            p *= 2
        masks = []
        if p == N:
            m = 1
            while m < N:
                masks.append(m)
                m <<= 1
        else:
            m = 1
            while m < p:
                masks.append(m)
                m <<= 1
            if p < N:                       # remainder fold step
                masks.append(p)
        total = 0.0
        for m in masks:
            mh = max(hop(r, r ^ m, k) for r in range(N))
            total += ser + mh * hop_lat
        return total

    if COLLECTIVE_ROUNDS.get(t) == "tree":
        # binomial tree depth, one serialized packet per level
        depth = max(1, int(math.ceil(math.log2(max(1, N)))))
        return depth * (ser + hop_lat * avg_hop(k))

    if t == "ALLTOALL":
        # each rank serially injects N-1 packets of bytes/N
        return (N - 1) * ser + hop_lat * avg_hop(k)

    return 0.0


def next_head_window(x, y, w, pw, t):
    """Earliest activation-window start >= t for mesh node (x, y) under the
    BookSim WSE phase clock (same geometry as TrafficManager::_wseAdvance /
    _wseUpdateStates):

      stage 0 (horizontal): phase p activates column {p, p+w, p+2w, ...}
                            -> node column x is HEAD at phase p = x % w
      stage 1 (vertical):   phase p activates row {p, p+w, p+2w, ...}
                            -> node row y is HEAD at phase p = y % w
    window start offset within period = (stage*w + p) * pw, period = 2*w*pw.

    Returns t unchanged when t already lies inside an activation window
    (the PE may then inject immediately, exactly like a HEAD router in the
    gated run); otherwise the earliest future window start. Monotone non-
    decreasing, so phased EST only ever pushes releases later (safe).
    """
    if pw <= 0 or w <= 1:
        return t  # degenerate: whole mesh HEAD, no gating effect
    period = 2 * w * pw
    r = t % period
    starts = [(x % w) * pw, w * pw + (y % w) * pw]
    for s in starts:
        if s <= r < s + pw:
            return t
    best = None
    for s in starts:
        cand = t + (s - r) % period
        if best is None or cand < best:
            best = cand
    return best


def node_dur(n, rate, freq_ratio, N, k, ser_lat, hop_lat):
    t = n.get("type", "")
    if t == "COMPUTE":
        if n.get("compute_ops", 0) > 0:
            return max(1.0, float(n["compute_ops"]) / rate)
        return max(1.0, float(n.get("compute_cycles", 0)) / freq_ratio)
    if t in COLLECTIVE_ROUNDS:
        return collective_dur(n, N, k, ser_lat, hop_lat)
    return 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ccdg", help="input CCDG JSON file")
    ap.add_argument("-o", "--out", default=None,
                    help="output prefix (default: input path without extension)")
    ap.add_argument("--cap", type=float, default=2.5e10,
                    help="compute capability ops/s (ccdg_compute_capability, "
                         "default 2.5e10 matching run_ccdg_mesh.sh)")
    ap.add_argument("--noc-ghz", type=float, default=2.0)
    ap.add_argument("--cpu-ghz", type=float, default=2.0)
    ap.add_argument("--flit-size", type=int, default=1,
                    help="flit_size_bytes (default 1, matching run_ccdg_mesh.sh)")
    ap.add_argument("--hop-lat", type=float, default=1.0,
                    help="per-hop latency in NoC cycles")
    ap.add_argument("--ser-lat", type=float, default=1.0,
                    help="per-flit serialization latency in NoC cycles")
    ap.add_argument("--gamma", type=float, default=1.0,
                    help="safety margin multiplier on flight/collective delays")
    ap.add_argument("--k", type=int, default=0,
                    help="mesh dimension (default int(sqrt(num_ranks)))")
    ap.add_argument("--pw", type=float, default=0.0,
                    help="wavelet phase width in NoC cycles; >0 enables "
                         "phased scheduling: every SEND/collective release is "
                         "aligned to the source rank's next HEAD window "
                         "(p=0: plain EST, responsive-equivalent)")
    ap.add_argument("--strip", type=int, default=0,
                    help="wavelet strip width w (0 = k). Same geometry as "
                         "wse_strip_width: phase p of the horizontal stage "
                         "activates columns {p,p+w,...}, w<=1 degenerates "
                         "(whole mesh HEAD, no effect)")
    args = ap.parse_args()

    d = json.load(open(args.ccdg))
    N = int(d["num_ranks"])
    nodes = d["nodes"]
    edges = d.get("cross_rank_edges", [])
    k = args.k or int(round(math.sqrt(N)))
    if k * k != N:
        print(f"WARNING: num_ranks={N} is not a perfect square; mesh k={k} "
              f"assumed (use --k to override)", file=sys.stderr)

    rate = args.cap / (args.noc_ghz * 1e9)
    freq_ratio = args.cpu_ghz / args.noc_ghz
    ser = args.ser_lat * args.gamma
    hlat = args.hop_lat * args.gamma

    # id -> index (ids are dense in practice, but stay safe)
    idx = {n["id"]: i for i, n in enumerate(nodes)}

    # per-node duration
    dur = [node_dur(n, rate, freq_ratio, N, k, ser, hlat) for n in nodes]
    M = len(nodes)

    # chain predecessors (same-rank, from 'predecessors' array)
    chain_prev = [-1] * M
    for i, n in enumerate(nodes):
        ps = n.get("predecessors", [])
        if not ps:
            continue
        # predecessors are chain-order; take the last one (the immediate one
        # on the rank's serial stream) and warrant est >= all of them
        chain_prev[i] = idx[ps[-1]]

    # cross edges: est[dst] >= est[src] + dur[src] + flight
    edge_delay = []
    for e in edges:
        s = idx[e["src_node"]]
        tgt = idx[e["dst_node"]]
        if s == tgt:
            continue
        bytes_ = nodes[s].get("comm_bytes", 0) or 64
        need = flits_of(bytes_, args.flit_size)
        flight = need * ser + hop(nodes[s]["rank"], nodes[tgt]["rank"], k) * hlat
        edge_delay.append((s, tgt, flight, need))

    # longest-path relaxation (DAG guarantees convergence)
    est = [0.0] * M
    changed = True
    iters = 0
    while changed and iters < 256:
        iters += 1
        changed = False
        # chain edges
        for i in range(M):
            p = chain_prev[i]
            if p < 0:
                continue
            v = est[p] + dur[p]
            if v > est[i]:
                est[i] = v
                changed = True
        # cross edges
        for s, tgt, flight, _ in edge_delay:
            v = est[s] + dur[s] + flight
            if v > est[tgt]:
                est[tgt] = v
                changed = True
    if iters >= 256:
        print("ERROR: relaxation did not converge (cyclic graph?)", file=sys.stderr)
        sys.exit(1)

    # ---- phased scheduling (wavelet-aligned EST) ----
    # The plain EST above is the responsive-equivalent lower bound; it never
    # staggers injections, so in-network queueing (blocked/congestion) is left
    # untouched. With --pw > 0 every SEND/collective release is pushed to the
    # source rank's next wavelet HEAD window, encoding the WSE phase clock
    # into the table: at release time the router is HEAD by construction, no
    # runtime gating is ever needed and contention windows are smoothed.
    # Fixed point: align -> re-propagate -> ... (est is monotone non-
    # decreasing and bounded, converges).
    w = args.strip if args.strip > 0 else k
    phased = args.pw > 0
    align_total = 0.0
    align_cnt = 0
    if phased:
        period = 2 * w * args.pw
        sys.stderr.write(
            f"[wavelet] phase width={args.pw:g} cyc strip={w} "
            f"period={period:g} cyc ({2 * w} windows)\n")
        for it in range(256):
            changed = False
            for i in range(M):
                t = nodes[i]["type"]
                if t not in ("SEND", "ISEND") and t not in COLLECTIVE_ROUNDS:
                    continue
                x = nodes[i]["rank"] % k
                y = nodes[i]["rank"] // k
                v = next_head_window(x, y, w, args.pw, est[i])
                if v > est[i]:
                    align_total += v - est[i]
                    align_cnt += 1
                    est[i] = v
                    changed = True
            if not changed:
                break
            for i in range(M):
                p = chain_prev[i]
                if p < 0:
                    continue
                v = est[p] + dur[p]
                if v > est[i]:
                    est[i] = v
            for s, tgt, flight, _ in edge_delay:
                v = est[s] + dur[s] + flight
                if v > est[tgt]:
                    est[tgt] = v

    crit = max(est[i] + dur[i] for i in range(M))

    # predicted orchestration slots: est gap along each rank's serial chain
    sched_slots = [0.0] * M
    total_slots = 0.0
    for i in range(M):
        p = chain_prev[i]
        if p < 0:
            sched_slots[i] = est[i]
        else:
            s = max(0.0, est[i] - (est[p] + dur[p]))
            sched_slots[i] = s
        total_slots += sched_slots[i]

    # ---- report ----
    types = collections.Counter(n["type"] for n in nodes)
    est_arr = sorted(est)
    print(f"== CCDG schedule report ==")
    print(f"ranks={N} mesh={k}x{k} nodes={M} cross_edges={len(edge_delay)}")
    print(f"compute rate={rate:.3f} ops/cyc (@cap {args.cap:g} ops/s, "
          f"noc {args.noc_ghz} GHz); ser={ser:g} cyc/flit hop_lat={hlat:g} cyc")
    print(f"node types: {dict(types)}")
    if phased:
        print(f"wavelet phased: pw={args.pw:g} cyc strip={w} "
              f"period={2 * w * args.pw:g} cyc; aligned {align_cnt} nodes, "
              f"total push = {align_total:,.0f} rank-cyc")
        est_arr2 = sorted(est)
        print(f"phased est distribution: min={est_arr2[0]:,.0f} "
              f"p50={est_arr2[len(est_arr2)//2]:,.0f} max={est_arr2[-1]:,.0f}")
    print(f"critical-path makespan bound = {crit:,.0f} cycles "
          f"({crit * 0.5e-3:.3f} ms @2GHz)")
    print(f"est distribution: min={est_arr[0]:,.0f} "
          f"p50={est_arr[len(est_arr)//2]:,.0f} max={est_arr[-1]:,.0f}")
    print(f"predicted orchestration slots (sched_wait bound) = "
          f"{total_slots:,.0f} rank-cyc "
          f"({total_slots / M:,.1f} /rank, {100.0 * total_slots / (M * crit):.1f}% of dwell)")
    if edge_delay:
        ds = sorted(x[2] for x in edge_delay)
        print(f"edge flight: n={len(ds)} min={ds[0]:,.0f} p50={ds[len(ds)//2]:,.0f} "
              f"max={ds[-1]:,.0f} (incl. {flits_of(max(x[2] for x in edge_delay), 1)} flit serialization)")

    # ---- outputs ----
    prefix = args.out if args.out else args.ccdg.rsplit(".", 1)[0]
    if phased:
        prefix = f"{prefix}_pw{args.pw:g}_w{w}"
    est_path = prefix + "_sched.est"
    with open(est_path, "w") as f:
        for i in range(M):
            f.write(f"{nodes[i]['id']} {est[i]:.0f}\n")
    print(f"wrote {est_path} ({M} lines)")

    # sorted CCDG: nodes reordered by (est, id), annotated with sched_est
    order = sorted(range(M), key=lambda i: (est[i], nodes[i]["id"]))
    out_nodes = []
    for i in order:
        n = dict(nodes[i])
        n["sched_est"] = round(est[i], 3)
        out_nodes.append(n)
    sorted_path = prefix + "_sorted.ccdg"
    with open(sorted_path, "w") as f:
        json.dump({"num_ranks": N, "nodes": out_nodes,
                   "cross_rank_edges": edges}, f)
    print(f"wrote {sorted_path}")


if __name__ == "__main__":
    main()