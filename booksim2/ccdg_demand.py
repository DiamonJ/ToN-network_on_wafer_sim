#!/usr/bin/env python3
"""Demand-driven mesh dataflow compiler (Cerebras paradigm).

Compile the LAMMPS COMMUNICATION DEMAND - the halo message set and the
compute phases between them - into a conflict-free mesh schedule. The MPI
trace is used only as the source of demand instances (who sends what to
whom, how much compute sits between phases); its timing and synchronization
structure (WAIT/WAITALL/WAITANY/RECV/IRECV/BARRIER) is DISCARDED.

Why: MPI synchronization is an artifact of the runtime, not of the physics.
Cerebras compiles the demand (stencil geometry, b+1 phases) so that data
arrives in a deterministic order and nothing ever waits. We do the same:
every SEND/collective injection reserves flit-exact slots on its DOR path
(LinkTable, destination ejection included), and happens-before constraints
raise each consumer's start to its data's arrival bound - so "data arrives
before it is consumed" is a compiler guarantee, not an assumption, and the
receiving rank needs no WAIT: blocked = 0 by construction. (--no-hb keeps
the legacy optimistic tier, which deletes the waits without satisfying
them; its makespan is a lower bound.)

Emitted artifacts:
  <prefix>_demand.ccdg  synchronization-free CCDG: per-rank chain of
                        COMPUTE/SEND/collective in trace order, no
                        cross_rank_edges, no WAIT*/RECV/BARRIER nodes
  <prefix>_demand.est   "node_id release_cycle" table consumed by BookSim's
                        ccdg_schedule_file channel (BookSim unchanged)

Model constants (identical to run_ccdg_mesh.sh / ccdg_scheduler):
  rate = cap / (noc_ghz * 1e9)        default 2.5e10 / 2e9 = 12.5 ops/cyc
  flit_size = 1, ser_lat = 1, hop_lat = 1
"""
import argparse
import collections
import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ccdg_scheduler import (  # noqa: E402
    COLLECTIVE_ROUNDS, flits_of, node_dur,
)
from ccdg_planner import LinkTable, dor_links  # noqa: E402


class PipelineLinkTable(LinkTable):
    """LinkTable variant with the BookSim IQ-router pipeline modeled:
    a flit occupies each hop's link for the transmission cycle, but the
    router pipeline (routing_delay=1 + vc_alloc_delay=1 + sw_alloc_delay=1
    + traverse=1) staggers consecutive hops by hop_stride cycles. Without
    this, the slot model underestimates link dwell and the real router
    queues up -> congestion. hop_stride=4 matches the default cfg."""

    def __init__(self, hop_stride=4.0):
        super().__init__()
        self.hop_stride = hop_stride
        self.retries = 0  # slot-conflict search restarts (assertion input)

    def first_free(self, links, flits, from_t, offsets=None):
        """Earliest injection t >= from_t such that every hop j occupies
        [t + offsets[j], t + offsets[j] + flits - 1] without overlap
        anywhere on its (pseudo-)link. offsets default to the plain
        hop-stride stagger; an appended ejection pseudo-link carries its
        own (larger) offset. Returns (t, arrival bound): the tail flit is
        fully at the destination at t + hops*stride + flits (ejection
        modeled: the ejection window ends at exactly that bound)."""
        hops = len(links)
        if offsets is None:
            offsets = [j * self.hop_stride for j in range(hops)]
        flight = hops * self.hop_stride + flits
        if len(offsets) > hops:  # ejection pseudo-link appended
            flight = offsets[-1] + flits
        cand = from_t
        self.conflict_owners = set()
        while True:
            esc = cand
            for lk, off in zip(links, offsets):
                a = cand + off
                b = a + flits - 1
                for (s, e, o) in self.intervals.get(lk, []):
                    if a <= e and b >= s:
                        esc = max(esc, e - off + 1)
                        self.conflict_owners.add(o)
            if esc == cand:
                return cand, flight
            cand = esc
            self.retries += 1
            if cand - from_t > (1 << 30):
                raise RuntimeError("first_free did not converge")

    def reserve(self, links, flits, t, offsets=None, owner=None):
        if offsets is None:
            offsets = [j * self.hop_stride for j in range(len(links))]
        for lk, off in zip(links, offsets):
            a = t + off
            b = t + off + flits - 1
            self.intervals[lk].append((a, b, owner))

# MPI runtime synchronization nodes: discarded entirely. RECV/IRECV only
# exist to set up later WAITs, so they carry no demand either.
SYNC_TYPES = frozenset(("WAIT", "WAITALL", "WAITANY", "RECV", "IRECV",
                        "BARRIER"))


# ----------------------------------------------------------------------
# module 1: demand extraction
# ----------------------------------------------------------------------

def extract_demand(path):
    """Walk a trace CCDG and keep only demand: per-rank chains of
    COMPUTE / SEND / collective, in trace order, without synchronization.

    Returns (N, k, chains, stats); chains[r] = [node dicts]; stats has
    conservation counters (comm_bytes by type, direction histogram).
    """
    d = json.load(open(path))
    N = int(d["num_ranks"])
    k = int(round(math.sqrt(N)))
    chains = [[] for _ in range(N)]
    dropped = collections.Counter()
    kept = collections.Counter()
    for nd in d["nodes"]:
        t = nd["type"]
        if t in SYNC_TYPES:
            dropped[t] += 1
            continue
        chains[nd["rank"]].append(nd)
        kept[t] += 1

    stats = {"dropped": dict(dropped), "kept": dict(kept),
             "comm_bytes_send": 0, "comm_bytes_coll": 0,
             "dirs": collections.Counter(),
             "dir_bytes": collections.Counter()}
    for r in range(N):
        for nd in chains[r]:
            if nd["type"] in ("SEND", "ISEND"):
                nb = nd.get("comm_bytes", 0) or 0
                stats["comm_bytes_send"] += nb
                s, t = nd["rank"], nd.get("comm_dst", (r + 1) % N)
                d = (t % k - s % k, t // k - s // k)
                stats["dirs"][d] += 1
                stats["dir_bytes"][d] += nb
            elif nd["type"] in COLLECTIVE_ROUNDS:
                stats["comm_bytes_coll"] += nd.get("comm_bytes", 0) or 0
    # Phase lattice: Cerebras "b controls which heads may inject".
    # b_d = DOR hop distance of direction d; w_d = b_d + 1 = footprint
    # width of one injected flow along d. Spacing same-direction flows
    # by w_d in the PERPENDICULAR axis separates their DOR paths in
    # space (disjoint link sets) -> zero contention independent of
    # flit count and timing.
    lattice = {}
    for d, cnt in stats["dirs"].items():
        bd = abs(d[0]) + abs(d[1])
        lattice[d] = {"b": bd, "w": bd + 1, "msgs": cnt,
                      "bytes": stats["dir_bytes"][d]}
    stats["lattice"] = lattice

    # Happens-before map: each SEND -> the first kept node on the receiving
    # rank after its matching WAIT (the consumer of that data). The demand
    # chain keeps trace order per rank, so the consumer is simply the next
    # kept node after the (dropped) WAIT in the receiving rank's sequence.
    # (d was rebound to a direction tuple by the lattice loop; reload.)
    data = json.load(open(path))
    ninfo = {nd["id"]: nd for nd in data["nodes"]}
    by_rank_all = collections.defaultdict(list)
    for nd in data["nodes"]:
        by_rank_all[nd["rank"]].append(nd["id"])
    kept_set = {nd["id"] for ch in chains for nd in ch}
    wait_sends = collections.defaultdict(list)
    for e in data.get("cross_rank_edges", []):
        a, b = e["src_node"], e["dst_node"]
        s, w = (a, b) if ninfo[a]["type"] in ("SEND", "ISEND") else (b, a)
        wait_sends[w].append(s)
    deps = []
    for w, sends in wait_sends.items():
        r = ninfo[w]["rank"]
        cid = next((nid for nid in by_rank_all[r]
                    if nid > w and nid in kept_set), None)
        if cid is not None:
            deps.extend((s, cid) for s in sends)
    stats["hb_deps"] = deps

    # Cerebras wavelet metadata per SEND: propagation stage (X if the first
    # DOR leg is horizontal, else Y), hop distance b_d, perpendicular-axis
    # coordinate c (y for horizontal moves, x for vertical moves - the
    # paper's strip index), and destination coordinates.
    send_meta = {}
    for r in range(N):
        for nd in chains[r]:
            if nd["type"] in ("SEND", "ISEND"):
                dst = nd.get("comm_dst", (r + 1) % N)
                dx = dst % k - r % k
                dy = dst // k - r // k
                b = abs(dx) + abs(dy)
                stage = "X" if dx != 0 else "Y"
                c = (r // k) if stage == "X" else (r % k)
                send_meta[nd["id"]] = {"stage": stage, "b": b, "c": c,
                                       "dst": dst, "dx": dx, "dy": dy}
    stats["send_meta"] = send_meta
    return N, k, chains, stats


# ----------------------------------------------------------------------
# module 1a: PBC fold renumbering (paper III-E interleaved placement)
# ----------------------------------------------------------------------

def fold_coord(i, k):
    """Split the periodic coordinate circle in two halves and collapse it
    to a line: first half ascending on even fabric coordinates, second
    half descending on odd ones. Ring neighbors land 1-2 hops apart
    (seams 1 hop, the rest 2), so no logical neighbor is ever more than
    2 hops away on the open mesh - the paper's "two hops away instead of
    one hop in the non-periodic case"."""
    return 2 * i if i < k // 2 else 2 * (k - 1 - i) + 1


def fold_rank(r, k):
    """Fold both dimensions independently: logical rank -> fabric node."""
    return fold_coord(r % k, k) + k * fold_coord(r // k, k)


def fold_renumber_ccdg(src, dst):
    """Renumber ranks (and communication endpoints) from the logical
    processor grid to the interleaved fabric layout. Node ids, cross-rank
    edges and est semantics are unchanged - only rank-space moves, so
    BookSim's rank==node invariant keeps holding with zero changes."""
    d = json.load(open(src))
    N = int(d["num_ranks"])
    k = int(round(math.sqrt(N)))

    def f(r):
        return fold_rank(r, k)

    for nd in d["nodes"]:
        if "rank" in nd:
            nd["rank"] = f(nd["rank"])
        for fld in ("comm_src", "comm_dst", "coll_root"):
            if fld in nd and isinstance(nd[fld], int):
                nd[fld] = f(nd[fld])
    with open(dst, "w") as fh:
        json.dump(d, fh)
    return k


# ----------------------------------------------------------------------
# module 1b: Cerebras stage-serialized chain reorder
# ----------------------------------------------------------------------

def reorder_chain(chains, N, send_meta):
    """Emit each rank's SENDs in Cerebras stage order: inside every
    COMPUTE/collective-delimited gap, X-stage messages come before
    Y-stage messages, each stage sorted by (perpendicular strip coord,
    destination) to match the wavefront sweep. COMPUTE anchors keep
    their chain positions (SENDs are fire-and-forget; consumers are
    guarded by HB constraints, not by chain adjacency)."""
    out = []
    for r in range(N):
        res, gap = [], []

        def flush():
            xs, ys = [], []
            for n in gap:
                (xs if send_meta[n["id"]]["stage"] == "X" else ys).append(n)
            key = lambda n: (send_meta[n["id"]]["c"],
                             send_meta[n["id"]]["dst"])
            xs.sort(key=key)
            ys.sort(key=key)
            res.extend(xs)
            res.extend(ys)

        for nd in chains[r]:
            if nd["type"] in ("SEND", "ISEND"):
                gap.append(nd)
            else:
                flush()
                gap = []
                res.append(nd)
        flush()
        out.append(res)
    return out


# ----------------------------------------------------------------------
# module 2: conflict-free slot planning (no cross-rank dependencies)
# ----------------------------------------------------------------------

def _rd_rounds(N):
    """Recursive-doubling round masks (identical to ccdg_scheduler's
    collective_dur and to CCDGTrafficManager::_collectivePattern)."""
    p = 1
    while p * 2 <= N:
        p *= 2
    if p == N:
        m = 1
        masks = []
        while m < N:
            masks.append(m)
            m <<= 1
    else:
        m = 1
        masks = []
        while m < p:
            masks.append(m)
            m <<= 1
        if p < N:
            masks.append(p)
    return masks


def plan_flow(chains, N, k, rate, freq_ratio, ser=1.0, hlat=1.0,
              slack=1.0, hop_stride=4.0, phase=False, lattice=None,
              phase_off=0, deps=None, ejection=True,
              cerebras=False, stage_barrier="global",
              wavelet_mode="serial", send_meta=None):
    """Slot every demand instance on the mesh with zero conflicts.

    Per rank two timelines:
      pe_t[r]  - PE chain progress (COMPUTE dwell, +1 per message node)
      port_t[r]- injection-port serialization end (packets leave one flit
                 per cycle; a later SEND may only inject after all earlier
                 packets have drained, otherwise BookSim backpressure would
                 stall it -> the wait we are eliminating).
    Release search is the earliest t >= max(pe_t, port_t) at which the whole
    DOR path is free (PipelineLinkTable, hop_stride = router pipeline).

    phase mode (Cerebras "b controls which heads may inject"): direction d
    messages may only inject at lattice points t ≡ c (mod w_d), where
    w_d = b_d + 1 and c is the PERPENDICULAR-axis coordinate of the source
    (x for horizontal moves, y for vertical moves - strips of width w_d,
    exactly the paper's horizontal/vertical stage strips). Same-direction
    flows on the same lattice are then separated in SPACE by >= w_d hops,
    so their DOR paths never share a link -> the slot search finds its
    first lattice point free on the first probe (zero retries, asserted).
    Off-lattice conflicts (other directions) are still absorbed by the
    slot table. The collective burst start snaps to the dominant w_d.

    Collective rounds inject back-to-back on the port (exactly like
    BookSim's non-blocking expansion), so they are placed as one
    contiguous burst - the whole burst shifts until every round's path is
    free. Round-robin single-node advance keeps the ranks fair.

    deps = [(send_id, consumer_id)]: happens-before constraints. A
    consumer node (the receiver's first demand node after its matching
    WAIT) may only be planned once every feeding SEND has a release time,
    and its start is raised to max(own chain, sender release + flight).
    This makes arrivals-before-consumption a compiler guarantee instead
    of an assumption. Ranks with unplanned feeders are deferred; the dep
    graph is a subgraph of the real execution's happens-before (acyclic).

    ejection: reserve the destination's ejection port as a pseudo-link at
    offset hops*stride, so simultaneous arrivals at one node serialize
    instead of silently overlapping (free-slot mode only; phase mode
    keeps its validated link-only assertion scope).

    Returns (plan_t {orig_id: cycle}, finish_t per rank,
    per_dir_retries Counter incl. 'collective' and 'total').
    """
    table = PipelineLinkTable(hop_stride=hop_stride)
    # cerebras pop-scan consumes the pending lists; work on a copy so the
    # caller's chains (emit, post-checks, reports) stay intact
    chains = [list(ch) for ch in chains]
    w_dom = 0
    if phase and lattice:
        w_dom = max(lp["w"] for lp in lattice.values())
    per_dir = collections.Counter()
    per_dir_same = collections.Counter()
    per_dir_cross = collections.Counter()
    pe_t = [0.0] * N
    port_t = [0.0] * N
    ptr = [0] * N
    plan_t = {}
    # happens-before: consumer node id -> feeding SEND node ids
    cons = collections.defaultdict(list)
    if deps:
        for s, c in deps:
            cons[c].append(s)
    node_by_id = {nd["id"]: nd for ch in chains for nd in ch}

    # Cerebras wavelet state: per-stage admission lattices with serial
    # rotation (phase q opens at the last landing of phase q-1) and an
    # X->Y stage barrier (global or per-rank). Y's phase 0 opens at the
    # barrier gate; X's phase 0 opens at 0.
    cb_tick = cb_remaining = cb_phase_open = None
    cb_x_left, cb_x_gmax = 0, 0.0
    cb_x_rank = cb_x_rank_max = None
    if cerebras and send_meta:
        ws = {"X": 1, "Y": 1}
        for m in send_meta.values():
            ws[m["stage"]] = max(ws[m["stage"]], m["b"] + 1)
        cb_stage = {sid: m["stage"] for sid, m in send_meta.items()}
        cb_tick = {sid: m["c"] % ws[cb_stage[sid]] for sid in send_meta}
        cb_remaining = {st: collections.Counter() for st in ws}
        cb_x_rank = [0] * N
        for ch in chains:
            for nd in ch:
                if nd["type"] in ("SEND", "ISEND") and nd["id"] in cb_stage:
                    st = cb_stage[nd["id"]]
                    cb_remaining[st][cb_tick[nd["id"]]] += 1
                    if st == "X":
                        cb_x_rank[nd["rank"]] += 1
        cb_x_left = sum(cb_remaining["X"].values())
        cb_phase_open = {st: [0.0] * ws[st] for st in ws}
        cb_x_rank_max = [0.0] * N

    def _stage_of(nd):
        if cerebras and send_meta and nd["id"] in send_meta:
            return cb_stage[nd["id"]], cb_tick[nd["id"]]
        return None, None

    def _ready(nd, r):
        """HB + Cerebras stage readiness for one node."""
        st, q = _stage_of(nd)
        if nd["id"] in cons and not all(s in plan_t for s in cons[nd["id"]]):
            return False, st, q
        if st is not None:
            if wavelet_mode == "serial" and q > 0 \
                    and cb_remaining[st][q - 1] > 0:
                return False, st, q
            if st == "Y":
                if stage_barrier == "global" and cb_x_left > 0:
                    return False, st, q
                if stage_barrier == "rank" and cb_x_rank[r] > 0:
                    return False, st, q
        if cerebras and st is None and stage_barrier != "none" \
                and nd["type"] in COLLECTIVE_ROUNDS:
            if stage_barrier == "global" and cb_x_left > 0:
                return False, st, q
            if stage_barrier == "rank" and cb_x_rank[r] > 0:
                return False, st, q
        return True, st, q

    any_active = True
    while any_active:
        any_active = False
        advanced = 0
        for r in range(N):
            if cerebras:
                # pop-scan: plan (and consume) the first ready node; a
                # barriered Y send cannot head-of-line block this rank's
                # later X sends, and a planned node is never re-planned
                pick, nd, nd_stage, nd_tick = -1, None, None, None
                for j in range(len(chains[r])):
                    ok, st_, q_ = _ready(chains[r][j], r)
                    if ok:
                        pick, nd, nd_stage, nd_tick = j, chains[r][j], st_, q_
                        break
                if pick < 0:
                    if chains[r]:
                        any_active = True  # deferred; retry next sweep
                    continue
                chains[r].pop(pick)
            else:
                if ptr[r] >= len(chains[r]):
                    continue
                nd = chains[r][ptr[r]]
                ok, nd_stage, nd_tick = _ready(nd, r)
                if not ok:
                    any_active = True
                    continue
                ptr[r] += 1
            any_active = True
            advanced += 1
            # arrival bound of the latest feeder: consumer must start at or
            # after the tail flit lands (t_inj + hops*stride + flits)
            hb_t = 0.0
            for s in cons.get(nd["id"], ()):
                snd = node_by_id[s]
                fl_s = flits_of(snd.get("comm_bytes", 0) or 64, 1)
                dst_s = snd.get("comm_dst", (snd["rank"] + 1) % N)
                hb_t = max(hb_t, plan_t[s]
                           + len(dor_links(snd["rank"], dst_s, k)[0]) * hop_stride
                           + fl_s)
            # Cerebras admission gate (final once deferral has passed)
            cb_gate = 0.0
            if nd_stage == "Y":
                if stage_barrier == "global":
                    cb_gate = cb_x_gmax
                elif stage_barrier == "rank":
                    cb_gate = cb_x_rank_max[r]
            elif nd_stage is None and cerebras and stage_barrier != "none" \
                    and nd["type"] in COLLECTIVE_ROUNDS:
                cb_gate = cb_x_gmax if stage_barrier == "global" \
                    else cb_x_rank_max[r]
            t = nd["type"]
            if t == "COMPUTE":
                if hb_t > pe_t[r]:
                    pe_t[r] = hb_t
                dur = node_dur(nd, rate, freq_ratio, N, k, ser, hlat)
                plan_t[nd["id"]] = pe_t[r]
                pe_t[r] += dur
            elif t in ("SEND", "ISEND"):
                dst = nd.get("comm_dst", (r + 1) % N)
                flits = flits_of(nd.get("comm_bytes", 0) or 64, 1)
                links, h_snd = dor_links(r, dst, k)
                t0 = max(pe_t[r], port_t[r], hb_t, cb_gate)
                ret_before = table.retries
                d = (dst % k - r % k, dst // k - r // k)
                if phase and lattice and d in lattice:
                    lp = lattice[d]
                    # admission lattice: t ≡ c (mod w_d); the phase axis is
                    # the axis perpendicular to the propagation axis
                    c = (r % k) if d[0] != 0 else (r // k)
                    c = (c + phase_off) % lp["w"]
                    w = lp["w"]
                    t0 = t0 + ((c - t0) % w)
                    while True:  # stay on the lattice across conflicts
                        t_inj, _ = table.first_free(links, flits,
                                                    from_t=t0)
                        if t_inj == t0:
                            break
                        # attribute the conflicts: same direction+phase
                        # violates the spatial-separation theorem, any
                        # other owner is normal lattice arbitration
                        same = sum(1 for o in table.conflict_owners
                                   if o and o[0] == d and o[1] == c)
                        per_dir_same[d] += same
                        per_dir_cross[d] += (len(table.conflict_owners)
                                             - same)
                        t0 = t_inj + ((c - t_inj) % w)
                    per_dir[d] += table.retries - ret_before
                    table.reserve(links, flits, t_inj, owner=(d, c))
                else:
                    if nd_stage is not None and wavelet_mode == "serial":
                        # serial rotation: this tick opened at the last
                        # landing of the previous tick (cb_gate carries the
                        # Y stage barrier)
                        t0 = max(t0, cb_phase_open[nd_stage][nd_tick])
                    elif nd_stage is not None and wavelet_mode == "interleave":
                        # control experiment: per-direction admission
                        # lattice, no serial rotation (barrier kept)
                        m = send_meta[nd["id"]]
                        wd = m["b"] + 1
                        c0 = (m["c"] + phase_off) % wd
                        t0 = t0 + ((c0 - t0) % wd)
                    # free-slot mode: model the destination ejection port
                    # as a pseudo-link at offset hops*stride so arrivals at
                    # the same node serialize instead of silently overlapping
                    offs = [j * hop_stride for j in range(len(links))]
                    e_links, e_offs = links, offs
                    if ejection:
                        e_links = links + [("ej", dst)]
                        e_offs = offs + [len(links) * hop_stride]
                    t_inj, _ = table.first_free(e_links, flits, from_t=t0,
                                                offsets=e_offs)
                    per_dir["free"] += table.retries - ret_before
                    table.reserve(e_links, flits, t_inj, offsets=e_offs,
                                  owner=("free", r))
                plan_t[nd["id"]] = t_inj
                if nd_stage is not None:
                    # wavelet bookkeeping: landing feeds the next tick's
                    # admission and the stage-barrier gates
                    land = t_inj + h_snd * hop_stride + flits
                    if nd_stage == "X":
                        cb_x_rank_max[r] = max(cb_x_rank_max[r], land)
                        cb_x_gmax = max(cb_x_gmax, land)
                        cb_x_left -= 1
                        cb_x_rank[r] -= 1
                    cb_remaining[nd_stage][nd_tick] -= 1
                    if nd_tick + 1 < len(cb_phase_open[nd_stage]):
                        cb_phase_open[nd_stage][nd_tick + 1] = max(
                            cb_phase_open[nd_stage][nd_tick + 1], land)
                pe_t[r] = max(pe_t[r], t_inj) + 1.0
                port_t[r] = t_inj + flits
            elif t in COLLECTIVE_ROUNDS:
                # BookSim expands all rounds at the arrival cycle and the
                # port drains them back-to-back; model the burst as one
                # contiguous reservation, shifting until every round path
                # is free.
                base = nd.get("comm_bytes", 0) or 64
                per = max(1, base // N)
                flits = flits_of(per, 1)
                rounds = _rd_rounds(N)
                t0 = max(pe_t[r], port_t[r], hb_t, cb_gate)
                ret_before = table.retries
                if phase and w_dom:
                    # burst start snaps to the dominant-direction lattice
                    c = ((r % k) + phase_off) % w_dom
                    t0 = t0 + ((c - t0) % w_dom)
                while True:  # find earliest contiguous burst start
                    cand = t0
                    ok = True
                    for m in rounds:
                        links, h_m = dor_links(r, r ^ m, k)
                        offs = [j * hop_stride for j in range(h_m)]
                        if ejection:
                            links = links + [("ej", r ^ m)]
                            offs = offs + [h_m * hop_stride]
                        # probe without reserving: whole path free from
                        # cand? (hop_stride dwell)
                        t_inj, _ = table.first_free(links, flits,
                                                    from_t=cand,
                                                    offsets=offs)
                        if t_inj > cand:
                            ok = False
                            t0 = t_inj
                            if phase and w_dom:
                                t0 = t0 + ((c - t0) % w_dom)
                            break
                        cand += flits
                    if ok:
                        break
                per_dir["collective"] += table.retries - ret_before
                cur = t0
                for m in rounds:
                    links, h_m = dor_links(r, r ^ m, k)
                    offs = [j * hop_stride for j in range(h_m)]
                    if ejection:
                        links = links + [("ej", r ^ m)]
                        offs = offs + [h_m * hop_stride]
                    table.reserve(links, flits, cur, offsets=offs,
                                  owner=("collective", r))
                    cur += flits
                plan_t[nd["id"]] = t0
                pe_t[r] = t0 + 1.0
                port_t[r] = cur
            else:  # unknown node type: zero-dwell passthrough
                plan_t[nd["id"]] = pe_t[r]
        if any_active and advanced == 0:
            raise RuntimeError("plan_flow deadlock: unsatisfiable "
                               "happens-before deps (acyclicity broken?)")
    finish = [pe_t[r] for r in range(N)]
    per_dir["total"] = table.retries
    return plan_t, finish, per_dir, per_dir_same, per_dir_cross


# ----------------------------------------------------------------------
# module 3: synchronization-free CCDG + release-time table
# ----------------------------------------------------------------------

def emit(chains, plan_t, N, out_prefix):
    """Renumber nodes, keep only demand fields, drop every dependency,
    annotate release times. cross_rank_edges is empty."""
    nodes = []
    for r in range(N):
        prev_id = None
        for nd in chains[r]:
            nid = len(nodes)
            n2 = {"id": nid, "rank": r, "type": nd["type"],
                  "sched_est": round(plan_t[nd["id"]], 3)}
            for f in ("compute_ops", "compute_cycles", "comm_src",
                      "comm_dst", "comm_bytes", "comm_count", "coll_root",
                      "wall_time_sec", "wall_duration_sec"):
                if f in nd:
                    n2[f] = nd[f]
            if prev_id is not None:
                n2["predecessors"] = [prev_id]
            nodes.append(n2)
            prev_id = nid
    out = {"num_ranks": N, "nodes": nodes, "cross_rank_edges": []}
    with open(out_prefix + "_demand.ccdg", "w") as f:
        json.dump(out, f)
    with open(out_prefix + "_demand.est", "w") as f:
        for nd in nodes:
            f.write("{} {}\n".format(nd["id"], nd["sched_est"]))
    return len(nodes)


# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ccdg", nargs="?", default=None,
                    help="input trace CCDG JSON (trimonly)")
    ap.add_argument("-o", "--out", default=None,
                    help="output prefix (default: input path w/o extension)")
    ap.add_argument("--fold-pbc", action="store_true",
                    help="paper III-E PBC fold: renumber ranks to the "
                         "interleaved fabric layout before planning, so "
                         "periodic-boundary neighbors are 1-2 hops away "
                         "instead of k-1 (requires even k)")
    ap.add_argument("--fold-only", nargs=2, metavar=("IN", "OUT"),
                    help="only apply the fold renumbering to a CCDG "
                         "(free-tier fold baseline) and exit")
    ap.add_argument("--cap", type=float, default=2.5e10,
                    help="compute capability ops/s (default 2.5e10)")
    ap.add_argument("--noc-ghz", type=float, default=2.0)
    ap.add_argument("--cpu-ghz", type=float, default=2.0)
    ap.add_argument("--hop-lat", type=float, default=1.0)
    ap.add_argument("--ser-lat", type=float, default=1.0)
    ap.add_argument("--slack", type=float, default=1.0,
                    help="flight-model margin (report only; release times "
                         "are slot-driven and need no slack)")
    ap.add_argument("--hop-stride", type=float, default=4.0,
                    help="router pipeline cycles per hop for link dwell "
                         "(routing+VA+SA+traverse = 4 in the default cfg)")
    ap.add_argument("--phase", action="store_true",
                    help="Cerebras b-parameterized injection admission: "
                         "direction d messages may only inject at lattice "
                         "points t ≡ c (mod w_d) with w_d = b_d + 1 "
                         "(default off = free-slot search)")
    ap.add_argument("--phase-offset", type=int, default=0,
                    help="lattice start offset added to the phase "
                         "coordinate (mod w_d)")
    ap.add_argument("--cerebras", action="store_true",
                    help="Cerebras wavelet compilation: reorder chains into "
                         "X/Y stages (destination-coordinate sweep order), "
                         "X->Y stage barrier, serial b+1 phase rotation; "
                         "happens-before constraints stay active")
    ap.add_argument("--stage-barrier", choices=["global", "rank", "none"],
                    default="global",
                    help="X->Y stage barrier scope (default global: the "
                         "whole grid finishes X before Y starts)")
    ap.add_argument("--wavelet-mode", choices=["serial", "interleave"],
                    default="serial",
                    help="serial = paper-style stage-internal phase "
                         "rotation (phase q opens at phase q-1's last "
                         "landing); interleave = per-direction admission "
                         "lattice, barrier kept (control experiment)")
    ap.add_argument("--no-hb", action="store_true",
                    help="disable happens-before release constraints "
                         "(legacy optimistic tier: consumers may start "
                         "before their data arrives - makespan is a "
                         "lower bound, not a correctness guarantee)")
    ap.add_argument("--no-ejection", action="store_true",
                    help="disable destination ejection-port modeling "
                         "(free-slot mode)")
    args = ap.parse_args()

    if args.fold_only:
        k = fold_renumber_ccdg(args.fold_only[0], args.fold_only[1])
        print(f"== fold-only: renumbered to interleaved fabric layout "
              f"(k={k}, paper III-E) -> {args.fold_only[1]}")
        return
    if not args.ccdg:
        ap.error("需要输入 CCDG 路径（或 --fold-only IN OUT）")

    prefix = args.out or (args.ccdg.rsplit(".", 1)[0] if "." in
                          args.ccdg.rsplit("/", 1)[-1]
                          else args.ccdg.rstrip("/") + "/demand")

    N, k, chains, stats = extract_demand(args.ccdg)
    if args.fold_pbc:
        if k % 2:
            print(f"ERROR: --fold-pbc 需要 k 为偶数 (k={k})", file=sys.stderr)
            sys.exit(2)
        fold_src = prefix + "_fold_input.ccdg"
        fold_renumber_ccdg(args.ccdg, fold_src)
        N, k, chains, stats = extract_demand(fold_src)
        bmax = max(m["b"] for m in stats["send_meta"].values()) \
            if stats["send_meta"] else 0
        n1 = sum(1 for m in stats["send_meta"].values() if m["b"] == 1)
        n2 = sum(1 for m in stats["send_meta"].values() if m["b"] == 2)
        print(f"== fold-pbc: logical->fabric rank 置换完成 (k={k}, "
              f"paper III-E): b1_msgs={n1} b2_msgs={n2} max_b={bmax}")
        if bmax > 2:
            print("ERROR: fold 后仍存在 b>2 方向（接缝未消解）",
                  file=sys.stderr)
            sys.exit(2)
    rate = args.cap / (args.noc_ghz * 1e9)
    freq_ratio = args.cpu_ghz / args.noc_ghz

    print("== demand extraction report ==")
    print(f"  num_ranks={N} mesh={k}x{k}")
    print(f"  kept nodes:   {dict(stats['kept'])}")
    print(f"  dropped sync: {stats['dropped']}")
    print(f"  comm bytes SEND={stats['comm_bytes_send']:,} "
          f"collective={stats['comm_bytes_coll']:,} "
          f"total={stats['comm_bytes_send'] + stats['comm_bytes_coll']:,}")
    print(f"  direction histogram (top 10): {stats['dirs'].most_common(10)}")
    if args.cerebras:
        chains = reorder_chain(chains, N, stats["send_meta"])
        hist = collections.Counter(
            (m["stage"], m["b"]) for m in stats["send_meta"].values())
        print("== cerebras stage serialization ==")
        print(f"  stage/b histogram: {dict(sorted(hist.items()))}")
    if args.phase:
        print("== phase lattice (Cerebras b -> injection-admission period) ==")
        for d, lp in sorted(stats["lattice"].items(),
                            key=lambda kv: -kv[1]["b"]):
            print(f"  dir {str(d):>9}  b={lp['b']:>2}  w={lp['w']:>2}  "
                  f"msgs={lp['msgs']:>5}  bytes={lp['bytes']:>12,}")
        w_dom = max(lp["w"] for lp in stats["lattice"].values())
        print(f"  dominant w_d = {w_dom} (collective burst-start lattice)")

    plan_t, finish, per_dir, per_dir_same, per_dir_cross = plan_flow(
        chains, N, k, rate, freq_ratio,
        ser=args.ser_lat, hlat=args.hop_lat,
        slack=args.slack,
        hop_stride=args.hop_stride,
        phase=args.phase, lattice=stats["lattice"],
        phase_off=args.phase_offset,
        deps=None if args.no_hb else stats.get("hb_deps"),
        ejection=(not args.no_ejection) and not args.phase,
        cerebras=args.cerebras, stage_barrier=args.stage_barrier,
        wavelet_mode=args.wavelet_mode,
        send_meta=stats.get("send_meta") if args.cerebras else None)

    # happens-before is a compiler guarantee: assert it post-plan
    hb_deps = stats.get("hb_deps") or []
    if hb_deps and not args.no_hb:
        node_by_id = {nd["id"]: nd for ch in chains for nd in ch}
        bad = 0
        for s, c in hb_deps:
            snd = node_by_id[s]
            fl_s = flits_of(snd.get("comm_bytes", 0) or 64, 1)
            dst_s = snd.get("comm_dst", (snd["rank"] + 1) % N)
            arr = (plan_t[s]
                   + len(dor_links(snd["rank"], dst_s, k)[0]) * args.hop_stride
                   + fl_s)
            if plan_t[c] < arr - 1e-6:
                bad += 1
        print(f"  happens-before: deps={len(hb_deps)} "
              f"consume-before-arrival violations={bad} (must be 0)")
        if bad:
            print("ERROR: happens-before violated after plan",
                  file=sys.stderr)
            sys.exit(2)
    tier = ("demand-legacy (no-hb, optimistic lower bound)" if args.no_hb
            else "demand+hb (correctness-guaranteed)")
    print(f"  tier: {tier}")
    if args.cerebras:
        sm = stats["send_meta"]
        node_by = {nd["id"]: nd for ch in chains for nd in ch}

        def _land(sid):
            m = sm[sid]
            snd = node_by[sid]
            return (plan_t[sid]
                    + len(dor_links(snd["rank"], m["dst"], k)[0]) * args.hop_stride
                    + flits_of(snd.get("comm_bytes", 0) or 64, 1))

        inv = 0
        ys = [sid for sid, m in sm.items() if m["stage"] == "Y"]
        if args.stage_barrier == "rank":
            xf = {}
            for sid, m in sm.items():
                if m["stage"] == "X":
                    r0 = node_by[sid]["rank"]
                    xf[r0] = max(xf.get(r0, 0.0), _land(sid))
            inv = sum(1 for sid in ys
                      if plan_t[sid] < xf.get(node_by[sid]["rank"], 0.0) - 1e-6)
        elif args.stage_barrier == "global":
            xf_all = max((_land(s) for s, m in sm.items()
                          if m["stage"] == "X"), default=0.0)
            inv = sum(1 for sid in ys if plan_t[sid] < xf_all - 1e-6)
        print(f"  stage barrier: {args.stage_barrier}, wavelet mode: "
              f"{args.wavelet_mode}, stage_inversions={inv} (must be 0)")
        if inv:
            print("ERROR: stage order violated after plan", file=sys.stderr)
            sys.exit(2)
        # pop-scan plans in bypass order; re-sort each rank's chain by
        # release time so the emitted chain order matches the est timeline
        # (BookSim walks the chain sequentially - a non-monotone est would
        # be unreachable and re-serialize the schedule)
        for r in range(N):
            chains[r].sort(key=lambda nd: plan_t[nd["id"]])
    makespan = max(finish)
    chain_sum = [sum(node_dur(nd, rate, freq_ratio, N, k, args.ser_lat,
                              args.hop_lat) for nd in chains[r])
                 for r in range(N)]
    print("== plan report ==")
    print(f"  makespan = {makespan:,.1f} cycles")
    print(f"  per-rank finish: min={min(finish):,.1f} "
          f"p50={sorted(finish)[N // 2]:,.1f} max={makespan:,.1f}")
    print(f"  per-rank pure chain length (zero-contention bound): "
          f"max={max(chain_sum):,.1f}")
    print(f"  contention slack absorbed = "
          f"{makespan - max(chain_sum):,.1f} cycles "
          f"({100.0 * (makespan - max(chain_sum)) / max(chain_sum):.1f}%)")
    if args.phase:
        print(f"  phase mode: on (offset={args.phase_offset}), "
              f"slot-conflict retries total={per_dir['total']}")
        print("  zero-contention check (same-direction, same-phase flows "
              "must be retry-free - spatial separation):")
        for d, lp in sorted(stats["lattice"].items(),
                            key=lambda kv: -kv[1]["b"]):
            same = per_dir_same.get(d, 0)
            cross = per_dir_cross.get(d, 0)
            axial = (d[0] == 0) != (d[1] == 0)
            if same == 0:
                status = "OK (spatially separated)"
            elif axial:
                status = "FAIL (axis-aligned dir must be retry-free)"
            else:
                status = "WARN (diagonal dir arbitrated by slot table)"
            print(f"    dir {str(d):>9}: same_phase_conflicts={same}  "
                  f"{status}; cross_phase_conflicts={cross} "
                  f"(slot-table arbitration)")
        print(f"    collective: retries={per_dir.get('collective', 0)}")

    m = emit(chains, plan_t, N, prefix)
    print(f"  emitted {m} nodes -> {prefix}_demand.ccdg + {prefix}_demand.est")


if __name__ == "__main__":
    main()
