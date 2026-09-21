#!/usr/bin/env python3
"""CCDG -> conflict-free deterministic plan (wait-free rendezvous orchestrator).

Upgrades the EST scheduler (ccdg_scheduler.py, zero-contention longest path)
to a RESOURCE-AWARE planner:

  * link time-slot table: every SEND/collective injection reserves the
    time ranges its flits occupy on each hop of its DOR path. Two messages
    never overlap on a link (no queuing -> congestion = 0 by construction).
  * release-time search: a message fires at the earliest cycle >= ready
    at which every hop of its path is free. Because the plan is conflict-free,
    every dependent data arrives BEFORE the consumer's scheduled slot ->
    blocked_wait -> 0 (the wait moves to the producer side: sched_wait).
  * balancing/shrink loop: the rendezvous point T' is a DECISION VARIABLE.
    Instead of anchoring on the slowest rank chain, bottleneck ranks get
    priority so fast ranks donate their slack; per-rank finish times are
    equalized and T' falls toward the est critical-path bound.
  * route rule: X-first-then-Y. Primary = DOR shortest path (identical to
    BookSim routing_function = dim_order). Detour candidates keep the same
    XY-first discipline: the X segment may end at an adjacent shim column,
    then the Y segment runs monotonically to the target row, and a tail X
    segment returns to the target column. NO YX-first path is ever emitted.
    (BookSim is untouched: physical packets still fly DOR, so a detour plan
    with longer flight - later release is CONSERVATIVE when injected.)
  * BARRIER: all ranks release their barrier node at the same planned
    cycle T_b (a rendezvous decision variable), not at the raw slowest
    arrival. Rebooking runs to a fixed point so barrier alignment never
    invalidates already reserved slots.

Model (identical constants to ccdg_scheduler / run_ccdg_mesh.sh):
    flit_size = 1, ser_lat = 1 cyc/flit, hop_lat = 1 cyc/hop
    arrival bound for a message = t_inj + hops + flits
Outputs:
    <prefix>_plan.est      "node_id plan_cycle" lines (consumed by BookSim
                           ccdg_schedule_file; per-node release-time gating)
    <prefix>_plan_sorted.ccdg
"""
import argparse
import collections
import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ccdg_scheduler import (  # noqa: E402
    COLLECTIVE_ROUNDS, flits_of, hop, next_head_window, node_dur,
)

# ----------------------------------------------------------------------
# link time-slot table
# ----------------------------------------------------------------------


def dor_links(src_rank, dst_rank, k):
    """DOR (X-first-then-Y) unidirectional links for rank src -> rank dst.

    Returns (links, hops); links are per-hop [(from(x,y), to(x,y))] tuples.
    """
    x1, y1 = src_rank % k, src_rank // k
    x2, y2 = dst_rank % k, dst_rank // k
    links = []
    x, y = x1, y1
    while x != x2:
        nx = x + (1 if x2 > x else -1)
        links.append(((x, y), (nx, y)))
        x = nx
    while y != y2:
        ny = y + (1 if y2 > y else -1)
        links.append(((x, y), (x, ny)))
        y = ny
    return links, len(links)


def xy_detour_links(src_rank, dst_rank, k, shim_col):
    """XY-first detour: X segment ends at column shim_col != x2, the Y
    segment runs monotonically to the target row, then a tail X segment
    returns to the target column (overall X-Y-X, never YX-first).

    Returns (links, hops) or None when degenerate / out of mesh.
    """
    x1, y1 = src_rank % k, src_rank // k
    x2, y2 = dst_rank % k, dst_rank // k
    if shim_col < 0 or shim_col >= k or shim_col == x2 or x1 == x2 or y1 == y2:
        return None
    links = []
    x, y = x1, y1
    while x != shim_col:
        nx = x + (1 if shim_col > x else -1)
        links.append(((x, y), (nx, y)))
        x = nx
    while y != y2:
        ny = y + (1 if y2 > y else -1)
        links.append(((x, y), (x, ny)))
        y = ny
    while x != x2:
        nx = x + (1 if x2 > x else -1)
        links.append(((x, y), (nx, y)))
        x = nx
    return links, len(links)


class LinkTable:
    """Per-unidirectional-link time interval registry (1 flit/cycle)."""

    def __init__(self):
        self.intervals = collections.defaultdict(list)  # link -> [(a, b)]

    def first_free(self, links, flits, from_t):
        """Earliest injection t >= from_t such that every hop j occupies the
        interval [t + j, t + j + flits - 1] without overlap anywhere.

        Returns (t, hops + flits): the arrival bound for the message."""
        hops = len(links)
        flight = hops + flits
        cand = from_t
        while True:
            esc = cand
            for (s, e) in self.intervals.get(links[0], []):
                a = cand
                b = cand + flits - 1
                if a <= e and b >= s:
                    esc = max(esc, e + 1)  # first free inject for hop 0
            if esc != cand:
                cand = esc
                continue
            for j in range(1, hops):
                lk = links[j]
                a = cand + j
                b = cand + j + flits - 1
                for (s, e) in self.intervals.get(lk, []):
                    if a <= e and b >= s:
                        esc = max(esc, e - j + 1)
                        break
                if esc > cand:
                    break
            if esc == cand:
                return cand, flight
            cand = esc
            if cand - from_t > (1 << 30):
                raise RuntimeError("first_free did not converge")

    def reserve(self, links, flits, t):
        for j, lk in enumerate(links):
            a = t + j
            b = t + j + flits - 1
            self.intervals[lk].append((a, b))


# ----------------------------------------------------------------------
# CCDG graph model
# ----------------------------------------------------------------------

class Ccdg:
    def __init__(self, path):
        d = json.load(open(path))
        self.n = int(d["num_ranks"])
        self.k = int(round(math.sqrt(self.n)))
        self.nodes = d["nodes"]
        self.m = len(self.nodes)
        self.idx = {nd["id"]: i for i, nd in enumerate(self.nodes)}
        self.edges = d.get("cross_rank_edges", [])

        self.chain_prev = [-1] * self.m
        for i, nd in enumerate(self.nodes):
            ps = nd.get("predecessors", [])
            if ps:
                self.chain_prev[i] = self.idx[ps[-1]]

        self.x_src = [self.idx[e["src_node"]] for e in self.edges]
        self.x_dst = [self.idx[e["dst_node"]] for e in self.edges]

    def dur(self, i, rate, freq_ratio, ser, hlat):
        return node_dur(self.nodes[i], rate, freq_ratio, self.n,
                        self.k, ser, hlat)


def base_est(g, rate, freq_ratio, ser, hlat, slack=1.0):
    """Zero-contention earliest-start (longest path on chain + cross edges).

    slack scales every flight time (conservative margin for BookSim
    router-pipeline delays the planner model does not capture).
    Returns (est, dur, flight_map) with flight_map {(s_idx, t_idx): cyc}."""
    m = g.m
    est = [0.0] * m
    dur = [g.dur(i, rate, freq_ratio, ser, hlat) for i in range(m)]
    flight_map = {}
    for s, t in zip(g.x_src, g.x_dst):
        if s == t:
            continue
        need = flits_of(g.nodes[s].get("comm_bytes", 0) or 64, 1)
        fly = (need * ser + hop(g.nodes[s]["rank"], g.nodes[t]["rank"],
                               g.k) * hlat) * slack
        flight_map[(s, t)] = fly
    for _ in range(256):
        changed = False
        for i in range(m):
            p = g.chain_prev[i]
            if p >= 0 and est[p] + dur[p] > est[i]:
                est[i] = est[p] + dur[p]
                changed = True
        for (s, t), fly in flight_map.items():
            v = est[s] + dur[s] + fly
            if v > est[t]:
                est[t] = v
                changed = True
        if not changed:
            break
    return est, dur, flight_map


# ----------------------------------------------------------------------
# the planner
# ----------------------------------------------------------------------

class Planner:
    def __init__(self, g, rate, freq_ratio, ser=1.0, hlat=1.0,
                 alt_paths=False, converge=8, rebook_max=8, slack=1.0):
        self.g = g
        self.rate = rate
        self.freq_ratio = freq_ratio
        self.ser = ser
        self.hlat = hlat
        self.alt_paths = alt_paths
        self.converge = max(1, converge)
        self.rebook_max = max(1, rebook_max)
        self.slack = slack

    # ---- message abstraction -------------------------------------
    def msg_units(self, i):
        """Expand node i into injection units [(src_rank, dst_rank, flits)].

        SEND/ISEND -> one unit. Collective -> one unit per expansion target
        (recursive doubling / binomial tree / alltoall), comm_bytes/N each.
        BARRIER -> [] (global rendezvous, handled separately).
        """
        nd = self.g.nodes[i]
        t = nd["type"]
        if t in ("SEND", "ISEND"):
            need = flits_of(nd.get("comm_bytes", 0) or 64, 1)
            return [(nd.get("comm_src", nd["rank"]),
                     nd.get("comm_dst", nd["rank"]), need)]
        if t in COLLECTIVE_ROUNDS:
            base = nd.get("comm_bytes", 0) or 64
            per = max(1, base // self.g.n)
            n = self.g.n
            r = nd["rank"]
            units = []
            if COLLECTIVE_ROUNDS[t] == "rd":
                p = 1
                while p * 2 <= n:
                    p *= 2
                masks = []
                msk = 1
                while msk < p:
                    masks.append(msk)
                    msk <<= 1
                if p < n:
                    masks.append(p)
                for mask in masks:
                    units.append((r, r ^ mask, per))
            elif COLLECTIVE_ROUNDS[t] == "tree":
                depth = max(1, int(math.ceil(math.log2(max(1, n)))))
                for lvl in range(depth):
                    d = r + (1 << lvl)
                    if t in ("BCAST", "SCATTER") and d < n:
                        units.append((r, d, per))
                    s = r - (1 << lvl)
                    if t in ("REDUCE", "GATHER") and s >= 0:
                        units.append((r, s, per))
            else:  # alltoall: serial fan-out
                for d in range(n):
                    if d != r:
                        units.append((r, d, per))
            return units
        return []

    # ---- single-message slot search ------------------------------
    def _seek(self, src, dst, flits, t0, table):
        """Earliest conflict-free injection. Primary DOR; with alt_paths the
        XY-first detours (shim columns x2-1 / x2+1) are also evaluated and
        the earliest feasible route wins (detour flight is longer, so the
        release is conservative for the real DOR flight).

        Returns (t, flight, links, detour) where links is the reserved path."""
        links, hops = dor_links(src, dst, self.g.k)
        t, flight = table.first_free(links, flits, t0)
        used, detour = links, False
        if self.alt_paths:
            x2 = dst % self.g.k
            for shim in (x2 - 1, x2 + 1):
                dl = xy_detour_links(src, dst, self.g.k, shim)
                if dl is None:
                    continue
                dlinks, dhops = dl
                t2, f2 = table.first_free(dlinks, flits, t0)
                if t2 < t:
                    t, flight, used, detour = t2, f2, dlinks, True
        return t, flight, used, detour

    # ---- one deterministic pass -----------------------------------
    def schedule_once(self, order, dur, flight_map, floor=None):
        """List-schedule nodes in order; returns (plan, book, finish, table).

        floor: previous plan (rebooking fixed point); every release is
        monotone non-decreasing across rebooks."""
        g = self.g
        m = g.m
        table = LinkTable()
        plan = [0.0] * m
        rank_busy = [0.0] * g.n
        book = []  # (t_inj, node_idx, src, dst, flits, arrival)

        for i in order:
            nd = g.nodes[i]
            r = nd["rank"]
            t_ready = floor[i] if floor is not None else 0.0
            p = g.chain_prev[i]
            if p >= 0:
                t_ready = max(t_ready, plan[p] + dur[p])
            for (s, fly) in self.in_flight.get(i, ()):
                t_ready = max(t_ready, plan[s] + dur[s] + fly)

            units = self.msg_units(i)
            if not units:
                plan[i] = t_ready
                continue

            min_t = max(t_ready, rank_busy[r])
            first = None
            for (src, dst, need) in units:
                ti, fly, links, det = self._seek(src, dst, need, min_t, table)
                if det:
                    self.n_detours += 1
                if first is None:
                    first = ti
                table.reserve(links, need, ti)
                rank_busy[r] = max(rank_busy[r], ti + need)
                min_t = max(min_t, ti + need)
                book.append((ti, i, src, dst, need, ti + fly))
            plan[i] = first

        # ---- BARRIER alignment: rendezvous point = decision variable
        b_idxs = [i for i in range(m) if g.nodes[i]["type"] == "BARRIER"]
        if b_idxs:
            t_b = max(plan[i] for i in b_idxs)
            for i in b_idxs:
                plan[i] = max(plan[i], t_b)

        # ---- monotone propagation (alignment may push tails later) --
        for _ in range(256):
            changed = False
            for i in range(m):
                p = g.chain_prev[i]
                if p >= 0 and plan[p] + dur[p] > plan[i]:
                    plan[i] = plan[p] + dur[p]
                    changed = True
            for tgt, inl in self.in_flight.items():
                for (s, fly) in inl:
                    v = plan[s] + dur[s] + fly
                    if v > plan[tgt]:
                        plan[tgt] = v
                        changed = True
            if not changed:
                break

        finish = [0.0] * g.n
        for i in range(m):
            r = g.nodes[i]["rank"]
            finish[r] = max(finish[r], plan[i] + dur[i])
        return plan, book, finish, table

    # ---- main -------------------------------------------------------
    def run(self):
        g, m = self.g, self.g.m
        self.n_detours = 0
        est, dur, flight_map = base_est(g, self.rate, self.freq_ratio,
                                        self.ser, self.hlat, self.slack)
        self.flight_map = flight_map
        # Index incoming edges by target node: schedule_once and the
        # propagation passes must not scan every edge per node (O(M*E)
        # explodes at 256r: ~311k nodes x ~47k edges).
        in_flight = {}
        for (s, t), fly in flight_map.items():
            in_flight.setdefault(t, []).append((s, fly))
        self.in_flight = in_flight

        # Initial pass: rank-major (chain-consistent) order, no floor.
        order0 = list(range(m))
        plan, book, finish, table = self.schedule_once(order0, dur, flight_map)

        # Bottleneck-first reorders (balancing): slowest rank gets priority;
        # each candidate order is rebooked to a fixed point so the plan is
        # conflict-free and monotone.
        best = (max(finish), max(finish) - min(finish), plan, book, finish, -1)
        floor = plan
        for it in range(1, self.converge):
            bottleneck = max(range(g.n), key=lambda r: finish[r])
            order = sorted(range(m),
                           key=lambda i: (0 if g.nodes[i]["rank"] == bottleneck
                                          else 1,
                                          -(est[i] + dur[i]), i))
            for _ in range(self.rebook_max):
                nb = self.schedule_once(order, dur, flight_map, floor=floor)
                nplan, nbook, nfinish, ntab = nb
                if all(nplan[i] == floor[i] for i in range(m)):
                    break
                floor = nplan
            plan, book, finish, table = nplan, nbook, nfinish, ntab
            mk, gp = max(finish), max(finish) - min(finish)
            if mk < best[0]:
                best = (mk, gp, plan, book, finish, bottleneck)

        self.makespan, self.gap, self.plan = best[0], best[1], best[2]
        self.book, self.finish = best[3], best[4]
        self.est = est
        self.dur = dur
        self.flight_map = flight_map
        return self


# ----------------------------------------------------------------------
# synthetic scenarios
# ----------------------------------------------------------------------

def synth_ccdg(kind, msg_flits=64, big_flits=512, long_cycles=76800,
               normal_cycles=9600, rounds=6):
    """Synthetic 16-rank (4x4) CCDG.

    kind='ring': zero-conflict wave (assertion 1) - every rank computes the
                 same amount in parallel, then sends two near-neighbour
                 messages whose shortest paths are disjoint; all ranks meet
                 at the trailing ALLREDUCE+BARRIER naturally, no waiting.
    kind='skew': bottleneck mismatch (assertion 2/3) - rank 0 computes 8x
                 longer and then fires one go message per rank; only after
                 receiving the go do the other ranks emit their two cross-
                 diameter messages, so all of them collide with rank 0's
                 512-flit message on the X corridor.

    Structure follows the 16r-long CCDG convention: each rank's node array
    holds its full chain and every cross-rank edge targets the receiver's
    WAITALL node.
    """
    if kind == "ring":
        return ring_ccdg(msg_flits, normal_cycles)
    return skew_ccdg(msg_flits, big_flits, long_cycles, normal_cycles,
                     rounds=rounds)


def ring_ccdg(msg_flits=64, normal_cycles=9600, n=16):
    nodes, edges = [], []
    nid = 0

    def add(rank, typ, pred=None, **kw):
        nonlocal nid
        nd = {"id": nid, "rank": rank, "type": typ}
        if pred is not None:
            nd["predecessors"] = [pred]
        nd.update(kw)
        nodes.append(nd)
        nid += 1
        return nid - 1

    s1_of, s2_of, wa_of = {}, {}, {}
    for r in range(n):
        c0 = add(r, "COMPUTE", compute_cycles=normal_cycles)
        s1 = add(r, "SEND", pred=c0, comm_src=r, comm_dst=(r + 1) % n,
                 comm_bytes=msg_flits, comm_count=msg_flits)
        s2 = add(r, "SEND", pred=s1, comm_src=r, comm_dst=(r + 3) % n,
                 comm_bytes=msg_flits // 4, comm_count=msg_flits // 4)
        # messages this rank receives: from (r-1) via SEND wave 1 and
        # from (r-3)%n via SEND wave 2
        i1 = add(r, "IRECV", pred=s2, comm_src=(r - 1) % n, comm_dst=r,
                 comm_bytes=msg_flits, comm_count=msg_flits)
        i2 = add(r, "IRECV", pred=i1, comm_src=(r - 3) % n, comm_dst=r,
                 comm_bytes=msg_flits // 4, comm_count=msg_flits // 4)
        wa = add(r, "WAITALL", pred=i2)
        wa_of[r] = wa
        s1_of[r], s2_of[r] = s1, s2
        ar = add(r, "ALLREDUCE", pred=wa, comm_bytes=128)
        add(r, "BARRIER", pred=ar)
    for r in range(n):
        edges.append({"src_node": s1_of[(r - 1) % n], "dst_node": wa_of[r]})
        edges.append({"src_node": s2_of[(r - 3) % n], "dst_node": wa_of[r]})
    return {"num_ranks": n, "nodes": nodes, "cross_rank_edges": edges}


def skew_ccdg(msg_flits=64, big_flits=512, long_cycles=76800,
               normal_cycles=9600, n=16, rounds=6):
    """Bottleneck-mismatch scenario.

    rank 0 computes long_cycles and then fires go SENDs (16 flits) to every
    other rank. After the go, all ranks run `rounds` exchange rounds: each
    round sends msg_flits (+ rank 0's big_flits) cross-diameter to
    (r+10)%n and msg_flits//4 to (r+2)%n, then waits for the two inbound
    messages - the queues build up round after round in free mode while a
    slot-planned schedule interleaves them deterministically.
    """
    nodes, edges = [], []
    nid = 0

    def add(rank, typ, pred=None, **kw):
        nonlocal nid
        nd = {"id": nid, "rank": rank, "type": typ}
        if pred is not None:
            nd["predecessors"] = [pred]
        nd.update(kw)
        nodes.append(nd)
        nid += 1
        return nid - 1

    go_of, wgo_of = {}, {}
    s1_of, s2_of, wa_of = [], [], []  # per (round, rank)
    for r in range(n):
        cyc = long_cycles if r == 0 else normal_cycles
        c0 = add(r, "COMPUTE", compute_cycles=cyc)
        if r == 0:
            # fire go SENDs to all ranks, then join the exchange rounds
            prev = c0
            for d in range(1, n):
                prev = add(0, "SEND", pred=prev, comm_src=0, comm_dst=d,
                           comm_bytes=16, comm_count=16)
                go_of[d] = prev
            chain = prev
        else:
            # wait for rank 0's go before the first exchange round
            ig = add(r, "IRECV", pred=c0, comm_src=0, comm_dst=r,
                     comm_bytes=16, comm_count=16)
            wg = add(r, "WAITALL", pred=ig)
            wgo_of[r] = wg
            chain = wg
        s1_of.append({})
        s2_of.append({})
        wa_of.append({})
        for rd in range(rounds):
            fl1 = big_flits if r == 0 else msg_flits
            s1 = add(r, "SEND", pred=chain, comm_src=r,
                     comm_dst=(r + 10) % n, comm_bytes=fl1,
                     comm_count=fl1)
            s2 = add(r, "SEND", pred=s1, comm_src=r,
                     comm_dst=(r + 2) % n, comm_bytes=fl1 // 4,
                     comm_count=fl1 // 4)
            i1 = add(r, "IRECV", pred=s2, comm_src=(r - 10) % n, comm_dst=r,
                     comm_bytes=msg_flits, comm_count=msg_flits)
            i2 = add(r, "IRECV", pred=i1, comm_src=(r - 2) % n, comm_dst=r,
                     comm_bytes=msg_flits // 4, comm_count=msg_flits // 4)
            wa = add(r, "WAITALL", pred=i2)
            s1_of[r][rd], s2_of[r][rd], wa_of[r][rd] = s1, s2, wa
            chain = wa
        ar = add(r, "ALLREDUCE", pred=chain, comm_bytes=128)
        add(r, "BARRIER", pred=ar)
    for r in range(1, n):
        edges.append({"src_node": go_of[r], "dst_node": wgo_of[r]})
    for rd in range(rounds):
        for r in range(n):
            edges.append({"src_node": s1_of[r][rd],
                          "dst_node": wa_of[(r + 10) % n][rd]})
            edges.append({"src_node": s2_of[r][rd],
                          "dst_node": wa_of[(r + 2) % n][rd]})
    return {"num_ranks": n, "nodes": nodes, "cross_rank_edges": edges}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("ccdg", nargs="?", help="input CCDG JSON (ignored with "
                                            "--synth)")
    ap.add_argument("-o", "--out", default=None,
                    help="output prefix (default: input without extension)")
    ap.add_argument("--synth", choices=["ring", "skew"], default=None,
                    help="generate synthetic scenario instead of reading CCDG")
    ap.add_argument("--cap", type=float, default=2.5e10)
    ap.add_argument("--noc-ghz", type=float, default=2.0)
    ap.add_argument("--cpu-ghz", type=float, default=2.5)
    ap.add_argument("--flit-size", type=int, default=1)
    ap.add_argument("--alt-paths", action="store_true",
                    help="evaluate XY-first detours in slot search")
    ap.add_argument("--converge", type=int, default=8)
    ap.add_argument("--msg-flits", type=int, default=64)
    ap.add_argument("--big-flits", type=int, default=512)
    ap.add_argument("--long-cycles", type=int, default=76800)
    ap.add_argument("--normal-cycles", type=int, default=9600)
    ap.add_argument("--rounds", type=int, default=6,
                    help="skew scenario: exchange rounds after rank 0's go")
    ap.add_argument("--slack", type=float, default=1.1,
                    help="flight-time relaxation factor (1.0 = tight model, "
                         ">1 covers BookSim router-pipeline overhead)")
    ap.add_argument("--hlat", type=float, default=4.0,
                    help="per-hop latency in the flight model (link=1 + "
                         "BookSim routing/vc/sw pipeline delays)")
    ap.add_argument("--pw", type=float, default=0.0,
                    help="wavelet phase width in NoC cycles; >0 overlays "
                         "phased gating: every SEND/collective release is "
                         "aligned to the source rank's next HEAD window "
                         "(0 = pure conflict-free plan)")
    ap.add_argument("--strip", type=int, default=0,
                    help="wavelet strip width w (0 = k); w<=1 degenerates "
                         "(whole mesh HEAD, no effect)")
    args = ap.parse_args()

    if args.synth:
        d = synth_ccdg(args.synth, args.msg_flits, args.big_flits,
                       args.long_cycles, args.normal_cycles, args.rounds)
        prefix = args.out or f"/tmp/synth_{args.synth}"
        ccdg_path = prefix + ".ccdg"
        with open(ccdg_path, "w") as f:
            json.dump(d, f)
        print(f"wrote synthetic CCDG: {ccdg_path} "
              f"({len(d['nodes'])} nodes, {len(d['cross_rank_edges'])} edges)")
        src = ccdg_path
    else:
        if not args.ccdg:
            ap.error("ccdg path or --synth is required")
        src = args.ccdg
        prefix = args.out or args.ccdg.rsplit(".", 1)[0]

    g = Ccdg(src)
    rate = args.cap / (args.noc_ghz * 1e9)
    freq_ratio = args.cpu_ghz / args.noc_ghz
    p = Planner(g, rate, freq_ratio, hlat=args.hlat, alt_paths=args.alt_paths,
                converge=args.converge, slack=args.slack).run()

    # ---- phased overlay (plan + wavelet): align every SEND/collective
    # release to the source rank's next HEAD window, then re-propagate the
    # chain + flight + BARRIER-alignment constraints to a fixed point.
    # All pushes are monotone non-decreasing, so the plan stays
    # conflict-free and the global rendezvous stays a decision variable.
    est_tag = ""
    if args.pw > 0:
        w = args.strip if args.strip > 0 else g.k
        est_tag = f"_pw{args.pw:g}_w{w}"
        period = 2 * w * args.pw
        plan = p.plan[:]
        align_total = 0.0
        align_cnt = 0
        b_idxs = [i for i in range(g.m)
                  if g.nodes[i]["type"] == "BARRIER"]
        for _ in range(256):
            changed = False
            for i in range(g.m):
                t = g.nodes[i]["type"]
                if t not in ("SEND", "ISEND") and t not in COLLECTIVE_ROUNDS:
                    continue
                x = g.nodes[i]["rank"] % g.k
                y = g.nodes[i]["rank"] // g.k
                v = next_head_window(x, y, w, args.pw, plan[i])
                if v > plan[i]:
                    align_total += v - plan[i]
                    align_cnt += 1
                    plan[i] = v
                    changed = True
            for i in range(g.m):
                pp = g.chain_prev[i]
                if pp >= 0 and plan[pp] + p.dur[pp] > plan[i]:
                    plan[i] = plan[pp] + p.dur[pp]
                    changed = True
            for (s, tgt), fly in p.flight_map.items():
                v = plan[s] + p.dur[s] + fly
                if v > plan[tgt]:
                    plan[tgt] = v
                    changed = True
            if b_idxs:
                t_b = max(plan[i] for i in b_idxs)
                for i in b_idxs:
                    if t_b > plan[i]:
                        plan[i] = t_b
                        changed = True
            if not changed:
                break
        p.plan = plan
        finish = [0.0] * g.n
        for i in range(g.m):
            r = g.nodes[i]["rank"]
            finish[r] = max(finish[r], plan[i] + p.dur[i])
        p.finish = finish
        p.makespan = max(finish)
        p.gap = max(finish) - min(finish)
        print(f"phased overlay: pw={args.pw:g} cyc strip={w} "
              f"period={period:g} cyc; aligned {align_cnt} nodes, "
              f"total push={align_total:,.0f} rank-cyc")

    types = collections.Counter(nd["type"] for nd in g.nodes)
    print("== CCDG conflict-free plan ==")
    print(f"ranks={g.n} mesh={g.k}x{g.k} nodes={g.m} "
          f"cross_edges={len(g.edges)} types={dict(types)}")
    est_cp = max(p.est[i] + p.dur[i] for i in range(g.m))
    print(f"zero-contention est CP = {est_cp:,.0f} cycles")
    print(f"detours={'on' if args.alt_paths else 'off'} "
          f"passes={args.converge} slack={args.slack:g} hlat={args.hlat:g}")
    print(f"makespan T' = {p.makespan:,.0f} cycles "
          f"({p.makespan * 0.5e-3:.3f} ms @2GHz)  "
          f"T'/estCP = {p.makespan / max(1.0, est_cp):.3f}")
    f0 = p.finish
    print(f"per-rank finish: min={min(f0):,.0f} max={max(f0):,.0f} "
          f"gap(max-min)={p.gap:,.0f} cycles "
          f"({p.gap / max(1.0, p.makespan) * 100:.2f}% of T')")
    if p.book:
        flights = sorted(b[5] - b[0] for b in p.book)
        print(f"injections={len(p.book)} flight p50="
              f"{flights[len(flights) // 2]:,.0f} max={flights[-1]:,.0f} "
              f"detours={p.n_detours}")

    est_path = prefix + est_tag + "_plan.est"
    with open(est_path, "w") as f:
        for i in range(g.m):
            f.write(f"{g.nodes[i]['id']} {p.plan[i]:.0f}\n")
    print(f"wrote {est_path} ({g.m} lines)")

    order = sorted(range(g.m), key=lambda i: (p.plan[i], g.nodes[i]["id"]))
    out_nodes = []
    for i in order:
        nd = dict(g.nodes[i])
        nd["sched_est"] = round(p.plan[i], 3)
        out_nodes.append(nd)
    spath = prefix + est_tag + "_plan_sorted.ccdg"
    with open(spath, "w") as f:
        json.dump({"num_ranks": g.n, "nodes": out_nodes,
                   "cross_rank_edges": g.edges}, f)
    print(f"wrote {spath}")


if __name__ == "__main__":
    main()