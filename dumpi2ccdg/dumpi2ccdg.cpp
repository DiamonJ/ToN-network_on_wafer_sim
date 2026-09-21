/**
 * dumpi2ccdg - Convert DUMPI trace files to Communication-Communication
 *              Dependency Graph (CCDG) in JSON format.
 *
 * Links against libundumpi for trace parsing.
 */
#include <dumpi/libundumpi/libundumpi.h>
#include <dumpi/common/argtypes.h>
#include <dumpi/common/constants.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cfloat>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <algorithm>
#include <glob.h>

/* ================================================================
 * Data structures
 * ================================================================ */

/* CPU-frequency assumption (ops/cycle basis) configurable via
 * CCDG_CPU_FREQ_GHZ. Legacy midpoint mode (CCDG_LEGACY_CPU_MID=1)
 * reproduces the old (start+stop)/2 extraction for regression against
 * compact_v2 (busy-wait self-spin was absorbed into COMPUTE there).
 * Declared before RankState because add_compute_node uses it. */
static double g_cpu_freq_ghz = 2.5;
static bool   g_legacy_cpu_mid = false;

enum NodeType { COMPUTE, COMM_SEND, COMM_RECV, COMM_ISEND, COMM_IRECV,
                COMM_WAIT, COMM_WAITANY, COMM_WAITALL, COMM_ALLREDUCE, COMM_BARRIER,
                COMM_BCAST, COMM_GATHER, COMM_ALLGATHER, COMM_SCATTER,
                COMM_ALLTOALL, COMM_REDUCE, COMM_OTHER };

static const char* node_type_str(NodeType t) {
    switch(t) {
        case COMPUTE:        return "COMPUTE";
        case COMM_SEND:      return "SEND";
        case COMM_RECV:      return "RECV";
        case COMM_ISEND:     return "ISEND";
        case COMM_IRECV:     return "IRECV";
        case COMM_WAIT:      return "WAIT";
        case COMM_WAITANY:   return "WAITANY";
        case COMM_WAITALL:   return "WAITALL";
        case COMM_ALLREDUCE: return "ALLREDUCE";
        case COMM_BARRIER:   return "BARRIER";
        case COMM_BCAST:     return "BCAST";
        case COMM_GATHER:    return "GATHER";
        case COMM_ALLGATHER: return "ALLGATHER";
        case COMM_SCATTER:   return "SCATTER";
        case COMM_ALLTOALL:  return "ALLTOALL";
        case COMM_REDUCE:    return "REDUCE";
        case COMM_OTHER:     return "OTHER";
        default:             return "UNKNOWN";
    }
}

struct CCDGNode {
    uint64_t    id;
    int         rank;
    NodeType    type;
    double      compute_cycles;   // CPU cycles spent in compute (between MPI calls)
    double      compute_time_sec; // CPU time in seconds
    double      compute_ops    = 0.0; // ops estimate (1 op/cycle @ CCDG_CPU_FREQ_GHZ)
    // Communication info (valid for COMM nodes)
    int         comm_src        = -1;
    int         comm_dst        = -1;
    int         comm_tag        = -1;
    uint64_t    comm_bytes      = 0;
    int         comm_count      = 0;
    int         comm_datatype_size = 0;
    // Collective info
    int         collective_root = -1;
    // Wall-clock timing (CLOCK_MONOTONIC seconds; for phase-trace correlation)
    double      wall_time_sec     = -1.0;  // midpoint of the MPI call
    double      wall_duration_sec = -1.0;  // stop - start of the MPI call
    // Dependency tracking
    std::vector<uint64_t> predecessors;
    // Non-blocking request tracking
    int         pending_req_id  = -1;   // -1 if not applicable
};

struct PendingReq {
    int         rank;
    int         req_id;          // local request identifier
    NodeType    type;            // ISEND or IRECV
    int         peer;            // source or destination rank
    int         tag;
    uint64_t    bytes;
    uint64_t    node_id;         // the CCDG node that created this request
};

/* ================================================================
 * Per-rank parsing state
 * ================================================================ */
struct RankState {
    int         rank_id;
    int         num_ranks;
    std::vector<CCDGNode> nodes;
    double      last_cpu_stop;   // last CPU time seen (seconds; gap end / legacy mid)
    uint64_t*   global_counter;  // pointer to global node counter
    int         next_req_id;
    std::map<int, PendingReq> pending_requests; // req_id -> PendingReq

    RankState() : rank_id(0), num_ranks(1), last_cpu_stop(0.0),
                  global_counter(0), next_req_id(0) {}

    /* Add a compute node between two MPI calls */
    void add_compute_node(double cpu_start, double cpu_stop) {
        double elapsed = cpu_stop - cpu_start;
        if (elapsed <= 0.0) return; // skip zero-length computes
        CCDGNode n;
        n.id = (*global_counter)++;
        n.rank = rank_id;
        n.type = COMPUTE;
        n.compute_time_sec = elapsed;
        n.compute_ops = elapsed * g_cpu_freq_ghz * 1e9; // ops estimate (1 op/cycle)
        n.compute_cycles = n.compute_ops;               // CPU cycles @ assumed freq
        n.comm_bytes = 0;
        n.pending_req_id = -1;
        // Predecessor is the previous node
        if (!nodes.empty())
            n.predecessors.push_back(nodes.back().id);
        nodes.push_back(n);
    }

    /* Add a communication node */
    uint64_t add_comm_node(NodeType type, double cpu_time,
                           int src=-1, int dst=-1, int tag=-1,
                           uint64_t bytes=0, int count=0,
                           int datatype_size=0, int req_id=-1) {
        CCDGNode n;
        n.id = (*global_counter)++;
        n.rank = rank_id;
        n.type = type;
        n.compute_time_sec = 0;
        n.compute_cycles = 0;
        n.comm_src = src;
        n.comm_dst = dst;
        n.comm_tag = tag;
        n.comm_bytes = bytes;
        n.comm_count = count;
        n.comm_datatype_size = datatype_size;
        n.pending_req_id = req_id;
        // Predecessor is the previous node
        if (!nodes.empty())
            n.predecessors.push_back(nodes.back().id);
        nodes.push_back(n);
        return n.id;
    }

    /* Add a collective node */
    uint64_t add_collective_node(NodeType type, double cpu_time,
                                 int bytes=0, int root=-1) {
        CCDGNode n;
        n.id = (*global_counter)++;
        n.rank = rank_id;
        n.type = type;
        n.compute_time_sec = 0;
        n.compute_cycles = 0;
        n.comm_bytes = bytes;
        n.collective_root = root;
        n.pending_req_id = -1;
        if (!nodes.empty())
            n.predecessors.push_back(nodes.back().id);
        nodes.push_back(n);
        return n.id;
    }

    /* Attach wall-clock timing to the node just added */
    void set_wall_time(double wt, double wd) {
        if (!nodes.empty()) {
            nodes.back().wall_time_sec = wt;
            nodes.back().wall_duration_sec = wd;
        }
    }
};

/* ================================================================
 * Global state across all ranks
 * ================================================================ */
struct GlobalState {
    int num_ranks;
    uint64_t global_node_counter;
    std::vector<RankState> ranks;
    // Cross-rank dependency edges
    std::vector<std::pair<uint64_t,uint64_t>> cross_edges; // (src_node_id, dst_node_id)

    GlobalState(int nr) : num_ranks(nr), global_node_counter(0), ranks(nr) {
        for (int i = 0; i < nr; i++) {
            ranks[i].rank_id = i;
            ranks[i].num_ranks = nr;
        }
    }
};

/* ================================================================
 * Callback helpers
 * ================================================================ */

static double cpu_start_to_sec(const dumpi_time* cpu) {
    return (double)cpu->start.sec + (double)cpu->start.nsec * 1e-9;
}

static double cpu_stop_to_sec(const dumpi_time* cpu) {
    return (double)cpu->stop.sec  + (double)cpu->stop.nsec  * 1e-9;
}

/* Legacy midpoint approximation of the call's CPU time */
static double cpu_mid_to_sec(const dumpi_time* cpu) {
    return (cpu_start_to_sec(cpu) + cpu_stop_to_sec(cpu)) / 2.0;
}

/* Emit the pure compute gap [last stop, this start) and return the
 * per-call CPU window (s, e) for the current MPI call. Legacy mode
 * keeps the old midpoint behavior: gap = [last mid, this mid). */
#define CPU_GAP_BEGIN(rs, cpu, s, e) \
    do { \
        if (g_legacy_cpu_mid) { \
            double _t = cpu_mid_to_sec(cpu); \
            (rs).add_compute_node((rs).last_cpu_stop, _t); \
            (s) = _t; (e) = _t; \
        } else { \
            (s) = cpu_start_to_sec(cpu); \
            (e) = cpu_stop_to_sec(cpu); \
            (rs).add_compute_node((rs).last_cpu_stop, (s)); \
        } \
    } while (0)

#define CPU_GAP_END(rs, e) do { (rs).last_cpu_stop = (e); } while (0)

static double wall_time_to_sec(const dumpi_time* wall) {
    double start = (double)wall->start.sec + (double)wall->start.nsec * 1e-9;
    double stop  = (double)wall->stop.sec  + (double)wall->stop.nsec  * 1e-9;
    return (start + stop) / 2.0;
}

static double wall_duration(const dumpi_time* wall) {
    double start = (double)wall->start.sec + (double)wall->start.nsec * 1e-9;
    double stop  = (double)wall->stop.sec  + (double)wall->stop.nsec  * 1e-9;
    return stop - start;
}

static double compute_datatype_size(dumpi_datatype dt) {
    // Basic MPI datatype sizes
    switch (dt) {
        case DUMPI_CHAR:           return 1;
        case DUMPI_SIGNED_CHAR:    return 1;
        case DUMPI_UNSIGNED_CHAR:  return 1;
        case DUMPI_BYTE:           return 1;
        case DUMPI_SHORT:          return 2;
        case DUMPI_UNSIGNED_SHORT: return 2;
        case DUMPI_INT:            return 4;
        case DUMPI_UNSIGNED:       return 4;
        case DUMPI_LONG:           return 8;
        case DUMPI_UNSIGNED_LONG:  return 8;
        case DUMPI_LONG_LONG:      return 8;
        case DUMPI_FLOAT:          return 4;
        case DUMPI_DOUBLE:         return 8;
        case DUMPI_LONG_DOUBLE:    return 16;
        default:                   return 4; // default to 4 bytes
    }
}

static int get_rank_from_userarg(void* userarg) {
    return static_cast<int>(reinterpret_cast<intptr_t>(userarg));
}

/* ================================================================
 * Generic callback (called for every MPI call)
 * ================================================================ */
struct CallbackContext {
    GlobalState* global;
    int rank;
};

/* ================================================================
 * Individual MPI call callbacks
 * ================================================================ */

static int on_send_cb(const dumpi_send* prm, uint16_t thread,
                      const dumpi_time* cpu, const dumpi_time* wall,
                      const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->count * compute_datatype_size(prm->datatype);
    uint64_t nid = rs.add_comm_node(COMM_SEND, s, ctx->rank, prm->dest,
                                    prm->tag, bytes, prm->count);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_recv_cb(const dumpi_recv* prm, uint16_t thread,
                      const dumpi_time* cpu, const dumpi_time* wall,
                      const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->count * compute_datatype_size(prm->datatype);
    uint64_t nid = rs.add_comm_node(COMM_RECV, s, prm->source, ctx->rank,
                                    prm->tag, bytes, prm->count);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_isend_cb(const dumpi_isend* prm, uint16_t thread,
                       const dumpi_time* cpu, const dumpi_time* wall,
                       const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->count * compute_datatype_size(prm->datatype);
    int req_id = rs.next_req_id++;
    uint64_t nid = rs.add_comm_node(COMM_ISEND, s, ctx->rank, prm->dest,
                                    prm->tag, bytes, prm->count, 0, req_id);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    PendingReq preq;
    preq.rank = ctx->rank;
    preq.req_id = req_id;
    preq.type = COMM_ISEND;
    preq.peer = prm->dest;
    preq.tag = prm->tag;
    preq.bytes = bytes;
    preq.node_id = nid;
    rs.pending_requests[req_id] = preq;
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_irecv_cb(const dumpi_irecv* prm, uint16_t thread,
                       const dumpi_time* cpu, const dumpi_time* wall,
                       const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->count * compute_datatype_size(prm->datatype);
    int req_id = rs.next_req_id++;
    uint64_t nid = rs.add_comm_node(COMM_IRECV, s, prm->source, ctx->rank,
                                    prm->tag, bytes, prm->count, 0, req_id);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    PendingReq preq;
    preq.rank = ctx->rank;
    preq.req_id = req_id;
    preq.type = COMM_IRECV;
    preq.peer = prm->source;
    preq.tag = prm->tag;
    preq.bytes = bytes;
    preq.node_id = nid;
    rs.pending_requests[req_id] = preq;
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_wait_cb(const dumpi_wait* prm, uint16_t thread,
                      const dumpi_time* cpu, const dumpi_time* wall,
                      const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    // Select the pending request this WAIT completes: prefer the oldest
    // IRECV (recv waits create cross-rank dependencies); fall back to the
    // oldest pending request of any type (send completion).
    int found_req_id = -1;
    for (auto& kv : rs.pending_requests) {
        if (kv.second.type == COMM_IRECV) { found_req_id = kv.first; break; }
    }
    if (found_req_id < 0 && !rs.pending_requests.empty()) {
        found_req_id = rs.pending_requests.begin()->first;
    }
    uint64_t nid = rs.add_comm_node(COMM_WAIT, s, -1, -1, -1, 0, 0, 0, found_req_id);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    // Add predecessor from the matching isend/irecv node
    if (found_req_id >= 0) {
        auto it = rs.pending_requests.find(found_req_id);
        if (it != rs.pending_requests.end()) {
            rs.nodes.back().predecessors.push_back(it->second.node_id);
            rs.pending_requests.erase(it);
        }
    }
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_waitall_cb(const dumpi_waitall* prm, uint16_t thread,
                         const dumpi_time* cpu, const dumpi_time* wall,
                         const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t nid = rs.add_comm_node(COMM_WAITALL, s);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    // Add predecessors from all pending isend/irecv nodes
    for (auto& kv : rs.pending_requests) {
        rs.nodes.back().predecessors.push_back(kv.second.node_id);
    }
    rs.pending_requests.clear();
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_waitany_cb(const dumpi_waitany* prm, uint16_t thread,
                         const dumpi_time* cpu, const dumpi_time* wall,
                         const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    // Same request-selection policy as MPI_Wait: prefer oldest IRECV
    int found_req_id = -1;
    for (auto& kv : rs.pending_requests) {
        if (kv.second.type == COMM_IRECV) { found_req_id = kv.first; break; }
    }
    if (found_req_id < 0 && !rs.pending_requests.empty()) {
        found_req_id = rs.pending_requests.begin()->first;
    }
    uint64_t nid = rs.add_comm_node(COMM_WAITANY, s, -1, -1, -1, 0, 0, 0, found_req_id);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    if (found_req_id >= 0) {
        auto it = rs.pending_requests.find(found_req_id);
        if (it != rs.pending_requests.end()) {
            rs.nodes.back().predecessors.push_back(it->second.node_id);
            rs.pending_requests.erase(it);
        }
    }
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_sendrecv_cb(const dumpi_sendrecv* prm, uint16_t thread,
                          const dumpi_time* cpu, const dumpi_time* wall,
                          const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    double wt = wall_time_to_sec(wall);
    double wd = wall_duration(wall);
    // Model as SEND + IRECV + WAIT so both directions participate in
    // cross-rank matching (the WAIT carries the recv dependency).
    uint64_t sbytes = prm->sendcount * compute_datatype_size(prm->sendtype);
    rs.add_comm_node(COMM_SEND, s, ctx->rank, prm->dest, prm->sendtag, sbytes, prm->sendcount);
    rs.set_wall_time(wt, wd);
    uint64_t rbytes = prm->recvcount * compute_datatype_size(prm->recvtype);
    int req_id = rs.next_req_id++;
    uint64_t rid = rs.add_comm_node(COMM_IRECV, s, prm->source, ctx->rank,
                                    prm->recvtag, rbytes, prm->recvcount, 0, req_id);
    rs.set_wall_time(wt, wd);
    PendingReq preq;
    preq.rank = ctx->rank;
    preq.req_id = req_id;
    preq.type = COMM_IRECV;
    preq.peer = prm->source;
    preq.tag = prm->recvtag;
    preq.bytes = rbytes;
    preq.node_id = rid;
    rs.pending_requests[req_id] = preq;
    uint64_t nid = rs.add_comm_node(COMM_WAIT, s, -1, -1, -1, 0, 0, 0, req_id);
    rs.set_wall_time(wt, wd);
    rs.nodes.back().predecessors.push_back(rid);
    rs.pending_requests.erase(req_id);
    CPU_GAP_END(rs, e);
    return 0;
}

/* ================================================================
 * Cross-rank edge construction (post pass over all parsed ranks)
 *
 * Chronological replay of all ranks' comm events by wall time:
 *   IRECV        -> posted to the rank's pending pool
 *   SEND / ISEND -> enqueued to the FIFO of key (src, dst, tag)
 *   WAIT / WAITANY / WAITALL -> completes one (WAIT) or all (WAITALL)
 *     pending IRECVs: prefer an IRECV whose matching send has already
 *     arrived (earliest send first) -- this mirrors real MPI completion
 *     order and keeps the dependency graph acyclic.
 * The completed IRECV is linked to the wait node (predecessor rewrite)
 * and the edge send_node -> wait_node is recorded.
 * ================================================================ */
static bool is_collective_type(NodeType t) {
    return t == COMM_ALLREDUCE || t == COMM_BARRIER || t == COMM_BCAST ||
           t == COMM_GATHER || t == COMM_ALLGATHER || t == COMM_SCATTER ||
           t == COMM_ALLTOALL || t == COMM_REDUCE;
}

/* Remove cross edges that participate in dependency cycles.
 *
 * The execution model treats every collective instance as a global
 * synchronization point (all ranks must reach instance i before any
 * rank can leave it). Under that model, a matched wait->send pair can
 * introduce a cycle when per-rank wall clocks are skewed, which would
 * deadlock the NoC simulator. We detect cycles over
 *   program-order chains + cross edges + same-instance collective links
 * and drop one cross edge per cycle (the wait then degrades to a
 * non-blocking wait, which is safe for liveness).
 */
static size_t eliminate_cycles(GlobalState& global) {
    // Compact node indexing: ids were assigned sequentially from 0
    size_t total_nodes = 0;
    for (int r = 0; r < global.num_ranks; r++)
        total_nodes += global.ranks[r].nodes.size();
    if (total_nodes == 0) return 0;

    // node_id -> rank and program position
    std::vector<int> node_rank(total_nodes, -1);
    std::vector<uint64_t> prog_next(total_nodes, UINT64_MAX);
    std::vector<int> coll_inst(total_nodes, -1);
    std::vector<std::vector<uint64_t>> inst_nodes;
    for (int r = 0; r < global.num_ranks; r++) {
        const RankState& rs = global.ranks[r];
        int inst = 0;
        for (size_t i = 0; i < rs.nodes.size(); i++) {
            uint64_t id = rs.nodes[i].id;
            if (id >= total_nodes) continue;
            node_rank[id] = r;
            if (i + 1 < rs.nodes.size()) prog_next[id] = rs.nodes[i + 1].id;
            if (is_collective_type(rs.nodes[i].type)) {
                if (inst >= (int)inst_nodes.size()) inst_nodes.emplace_back();
                coll_inst[id] = inst;
                inst_nodes[inst].push_back(id);
                inst++;
            }
        }
    }

    size_t removed_total = 0;
    // Contract each collective instance (one node per rank) into a single
    // supernode: entering a collective is a global sync point, so all its
    // per-rank nodes are one logical node. Supernode ids start at total_nodes.
    size_t n_super = inst_nodes.size();
    size_t ext_nodes = total_nodes + n_super;
    // map any node id to its graph id (collective nodes -> supernode)
    std::vector<uint64_t> gmap(total_nodes);
    for (uint64_t u = 0; u < total_nodes; u++) gmap[u] = u;
    for (size_t k = 0; k < n_super; k++)
        for (uint64_t m : inst_nodes[k]) gmap[m] = total_nodes + k;

    for (int iter = 0; iter < 10000; iter++) {
        // Build adjacency: program chain + current cross edges (collective
        // instances contracted into supernodes)
        std::vector<std::vector<uint64_t>> adj(ext_nodes);
        for (uint64_t u = 0; u < total_nodes; u++) {
            if (prog_next[u] != UINT64_MAX)
                adj[gmap[u]].push_back(gmap[prog_next[u]]);
        }
        // From a supernode, continue on every rank that has this instance
        for (size_t k = 0; k < n_super; k++) {
            for (uint64_t m : inst_nodes[k]) {
                if (prog_next[m] != UINT64_MAX)
                    adj[total_nodes + k].push_back(gmap[prog_next[m]]);
            }
        }
        for (const auto& e : global.cross_edges) {
            if (e.first < total_nodes && e.second < total_nodes)
                adj[gmap[e.first]].push_back(gmap[e.second]);
        }

        // Iterative DFS to find one cycle
        std::vector<char> color(ext_nodes, 0); // 0 white 1 gray 2 black
        std::vector<uint64_t> parent(ext_nodes, UINT64_MAX);
        uint64_t cyc_u = UINT64_MAX, cyc_v = UINT64_MAX; // back edge u -> v
        for (uint64_t start = 0; start < ext_nodes && cyc_u == UINT64_MAX; start++) {
            if (color[start]) continue;
            std::vector<std::pair<uint64_t, size_t>> stack;
            color[start] = 1;
            stack.push_back({start, 0});
            while (!stack.empty() && cyc_u == UINT64_MAX) {
                uint64_t u = stack.back().first;
                size_t& idx = stack.back().second;
                if (idx < adj[u].size()) {
                    uint64_t v = adj[u][idx++];
                    if (v == u) continue; // contracted self-loop
                    if (color[v] == 1) { cyc_u = u; cyc_v = v; break; }
                    if (color[v] == 0) {
                        color[v] = 1; parent[v] = u;
                        stack.push_back({v, 0});
                    }
                } else {
                    color[u] = 2;
                    stack.pop_back();
                }
            }
        }
        if (cyc_u == UINT64_MAX) break; // acyclic

        // Reconstruct the cycle via parent chain (graph ids) and find a
        // cross edge whose contracted endpoints lie on it
        std::set<uint64_t> on_cycle;
        uint64_t x = cyc_u;
        on_cycle.insert(x);
        while (x != cyc_v && parent[x] != UINT64_MAX) {
            x = parent[x];
            on_cycle.insert(x);
        }
        int drop_idx = -1;
        for (size_t i = 0; i < global.cross_edges.size(); i++) {
            const auto& e = global.cross_edges[i];
            if (e.first >= total_nodes || e.second >= total_nodes) continue;
            if (on_cycle.count(gmap[e.first]) && on_cycle.count(gmap[e.second])) {
                drop_idx = (int)i;
                break;
            }
        }
        if (drop_idx < 0) {
            // Cycle without a cross edge (should not happen); stop to avoid looping
            fprintf(stderr, "WARNING: cycle without cross edge; stopping cycle elimination\n");
            break;
        }
        global.cross_edges.erase(global.cross_edges.begin() + drop_idx);
        removed_total++;
    }
    if (removed_total > 0) {
        fprintf(stderr, "Acyclicity pass: removed %zu cross edges that formed cycles under collective-sync semantics\n", removed_total);
    }
    return removed_total;
}

static void build_cross_edges(GlobalState& global) {
    struct Ev {
        double wall;
        int rank;
        int kind;       // 0=irecv 1=send 2=wait 3=waitall
        uint64_t node_id;
        int idx;        // index in rank's node vector
        std::tuple<int,int,int> key;
    };
    std::vector<Ev> events;

    for (int r = 0; r < global.num_ranks; r++) {
        const RankState& rs = global.ranks[r];
        for (size_t i = 0; i < rs.nodes.size(); i++) {
            const CCDGNode& n = rs.nodes[i];
            if (n.type == COMM_IRECV) {
                events.push_back({n.wall_time_sec, r, 0, n.id, (int)i,
                                  std::make_tuple(n.comm_src, r, n.comm_tag)});
            } else if (n.type == COMM_SEND || n.type == COMM_ISEND) {
                events.push_back({n.wall_time_sec, r, 1, n.id, (int)i,
                                  std::make_tuple(r, n.comm_dst, n.comm_tag)});
            } else if (n.type == COMM_WAIT || n.type == COMM_WAITANY) {
                events.push_back({n.wall_time_sec, r, 2, n.id, (int)i, {}});
            } else if (n.type == COMM_WAITALL) {
                events.push_back({n.wall_time_sec, r, 3, n.id, (int)i, {}});
            }
        }
    }
    std::sort(events.begin(), events.end(),
              [](const Ev& a, const Ev& b) { return a.wall < b.wall; });

    struct Posted { std::tuple<int,int,int> key; uint64_t irecv_id; int irecv_idx; };
    std::vector<std::vector<Posted>> pending(global.num_ranks); // outstanding IRECVs
    std::map<std::tuple<int,int,int>, std::deque<uint64_t>> send_q; // arrived sends

    size_t n_edges = 0, n_unmatched = 0;
    for (const Ev& ev : events) {
        if (ev.kind == 0) {
            pending[ev.rank].push_back({ev.key, ev.node_id, ev.idx});
        } else if (ev.kind == 1) {
            send_q[ev.key].push_back(ev.node_id);
        } else {
            // wait / waitany / waitall: complete pending IRECVs of this rank
            std::vector<size_t> to_complete;
            if (ev.kind == 3) {
                for (size_t i = 0; i < pending[ev.rank].size(); i++) to_complete.push_back(i);
            } else {
                // choose the pending IRECV whose send arrived earliest
                int best = -1; double best_wall = 0.0;
                for (size_t i = 0; i < pending[ev.rank].size(); i++) {
                    auto qit = send_q.find(pending[ev.rank][i].key);
                    if (qit != send_q.end() && !qit->second.empty()) {
                        // head-of-queue send wall time unknown here; use posting
                        // order among ready candidates (FIFO per key is enough)
                        if (best < 0) { best = (int)i; }
                    }
                }
                (void)best_wall;
                if (best < 0 && !pending[ev.rank].empty()) best = 0; // send not seen yet
                if (best >= 0) to_complete.push_back((size_t)best);
            }
            // complete in reverse index order so erase() keeps indices valid
            std::sort(to_complete.begin(), to_complete.end(), std::greater<size_t>());
            for (size_t ci : to_complete) {
                Posted p = pending[ev.rank][ci];
                pending[ev.rank].erase(pending[ev.rank].begin() + ci);
                auto qit = send_q.find(p.key);
                if (qit != send_q.end() && !qit->second.empty()) {
                    global.cross_edges.push_back({qit->second.front(), ev.node_id});
                    qit->second.pop_front();
                    n_edges++;
                } else {
                    n_unmatched++;
                }
                // rewrite the wait node's predecessor to the completed IRECV
                CCDGNode& wnode = global.ranks[ev.rank].nodes[ev.idx];
                CCDGNode& inode = global.ranks[ev.rank].nodes[p.irecv_idx];
                bool has = false;
                for (uint64_t pid : wnode.predecessors) if (pid == p.irecv_id) has = true;
                if (!has) wnode.predecessors.push_back(p.irecv_id);
                (void)inode;
            }
        }
    }

    size_t leftover_sends = 0;
    for (auto& kv : send_q) leftover_sends += kv.second.size();
    size_t leftover_recvs = 0;
    for (auto& pv : pending) leftover_recvs += pv.size();
    fprintf(stderr, "Cross-edge matching (chronological): edges=%zu unmatched_waits=%zu leftover_sends=%zu leftover_posted_recvs=%zu\n",
            n_edges, n_unmatched, leftover_sends, leftover_recvs);
}

/* ================================================================
 * Setup trimming & burst compaction (optional post-processing)
 * ================================================================ */

static bool is_p2p_type(NodeType t) {
    return t == COMM_SEND || t == COMM_ISEND || t == COMM_RECV || t == COMM_IRECV ||
           t == COMM_WAIT || t == COMM_WAITANY || t == COMM_WAITALL;
}

static bool is_burst_type(NodeType t) {
    return t == COMPUTE || is_p2p_type(t);
}

/* Renumber all node ids to be contiguous from 0; update predecessors
 * and cross edges accordingly. */
static void renumber_ids(GlobalState& global) {
    std::map<uint64_t, uint64_t> idmap;
    uint64_t next = 0;
    for (auto& rs : global.ranks)
        for (auto& n : rs.nodes) idmap[n.id] = next++;
    for (auto& rs : global.ranks) {
        for (auto& n : rs.nodes) {
            n.id = idmap[n.id];
            for (auto& p : n.predecessors) p = idmap[p];
        }
    }
    for (auto& e : global.cross_edges) {
        e.first = idmap[e.first];
        e.second = idmap[e.second];
    }
    global.global_node_counter = next;
}

/* Detect the setup/run boundary of one rank's node sequence.
 * Primary: the first P2P node, if the prefix is purely COMPUTE +
 * collectives (LAMMPS setup broadcasts; halo exchange starts in run).
 * Fallback: steady-state period search on the type sequence — the tail
 * window must contain enough P2P nodes and agree at period p with the
 * window one step earlier; boundary = earliest matching window start,
 * keeping at least 3 steady periods.
 * Returns 0 when undecided (no trim). */
static size_t detect_setup_boundary(const std::vector<CCDGNode>& nodes) {
    const size_t len = nodes.size();
    if (len < 16) return 0;
    std::vector<int> seq(len);
    for (size_t i = 0; i < len; i++) seq[i] = (int)nodes[i].type;

    // Primary: setup prefix = COMPUTE + collectives only (no P2P)
    {
        size_t first_p2p = len;
        for (size_t i = 0; i < len; i++)
            if (is_p2p_type(nodes[i].type)) { first_p2p = i; break; }
        if (first_p2p > 0 && first_p2p < len) {
            size_t coll = 0, comp = 0;
            bool clean = true;
            for (size_t i = 0; i < first_p2p; i++) {
                if (is_collective_type(nodes[i].type)) coll++;
                else if (nodes[i].type == COMPUTE) comp++;
                else { clean = false; break; }
            }
            if (clean && coll >= 1 && comp >= 1) return first_p2p;
        }
    }

    // Fallback: steady-state period search
    for (size_t p = 8; p <= 2000 && 8 * p <= len; p++) {
        // steady window = last 2 periods; require P2P density of an MD loop
        size_t wstart = len - 2 * p;
        size_t p2p_cnt = 0;
        for (size_t i = wstart; i < len; i++)
            if (is_p2p_type(nodes[i].type)) p2p_cnt++;
        if (p2p_cnt < p / 4) continue;
        size_t tot = 0, ok = 0;
        for (size_t i = wstart; i + p < len; i++) {
            tot++;
            if (seq[i] == seq[i + p]) ok++;
        }
        if (tot == 0 || (double)ok / (double)tot <= 0.95) continue;
        // boundary = earliest start of a window matching the steady window;
        // always keep at least 3 steady periods after trimming
        for (size_t s = 0; s + 3 * p <= len; s++) {
            size_t t2 = 0, ok2 = 0;
            for (size_t q = 0; q < p; q++) {
                t2++;
                if (seq[s + q] == seq[wstart + q]) ok2++;
            }
            if (t2 > 0 && (double)ok2 / (double)t2 > 0.9) return s;
        }
    }
    return 0;
}

/* Parse the total "Loop time of X on N procs for M steps" from the
 * LAMMPS log in trace_dir (log.lammps or lammps.log). Returns -1 when
 * unavailable. Used to anchor the steady-state window: one-shot costs
 * of a captured run (first neighbor build, init exchanges AFTER the
 * first P2P, timing printout tail) live outside this window and must
 * be trimmed, otherwise per-rank compute sums inflate several-fold and
 * SimGrid frequency calibration diverges beyond physical limits. */
static double read_loop_time(const std::string& trace_dir) {
    const char* names[] = {"log.lammps", "lammps.log"};
    for (const char* name : names) {
        std::ifstream f(trace_dir + "/" + name);
        if (!f.is_open()) continue;
        std::string line;
        double lt = -1.0;
        while (std::getline(f, line)) {
            size_t p = line.find("Loop time of");
            if (p == std::string::npos) continue;
            lt = atof(line.c_str() + p + 12); // keep the last occurrence
        }
        if (lt > 0.0) return lt;
    }
    return -1.0;
}

/* Trim the setup prefix of every rank (env CCDG_TRIM_SETUP).
 * Two stages per rank:
 *   1. detect_setup_boundary: drop the pre-first-P2P setup prefix.
 *   2. Step window (when the LAMMPS log is available): keep only the
 *      single captured step. LAMMPS Timer sync mode brackets the run
 *      loop with MPI_Barrier (barrier_start / barrier_stop), which are
 *      exactly the endpoints of "Loop time" in log.lammps; the step
 *      window is therefore [BARRIER-before-last-P2P,
 *      BARRIER-after-last-P2P]. The setup prefix (incl. neighbor build
 *      after the first P2P) and the post-step thermo printout
 *      (ALLREDUCE burst + BCAST + print COMPUTE) both lie outside this
 *      pair and are trimmed. Falls back to a Loop-time wall window
 *      when no such BARRIER pair exists. */
static void trim_setup(GlobalState& global, const std::string& trace_dir) {
    size_t trimmed_total = 0;
    for (int r = 0; r < global.num_ranks; r++) {
        RankState& rs = global.ranks[r];
        size_t b = detect_setup_boundary(rs.nodes);
        if (b == 0 || b >= rs.nodes.size()) {
            fprintf(stderr, "Trim setup: rank %d no boundary detected, keeping all %zu nodes\n",
                    r, rs.nodes.size());
            continue;
        }
        rs.nodes.erase(rs.nodes.begin(), rs.nodes.begin() + b);
        rs.nodes.front().predecessors.clear();
        trimmed_total += b;
        fprintf(stderr, "Trim setup: rank %d removed %zu setup nodes (%zu remain)\n",
                r, b, rs.nodes.size());
    }

    // Stage 2: Loop-time window trim
    const double loop_time = read_loop_time(trace_dir);
    if (loop_time > 0.0) {
        fprintf(stderr, "Trim setup: loop-time window enabled (Loop time = %.6f s)\n", loop_time);
        size_t window_total = 0;
        for (int r = 0; r < global.num_ranks; r++) {
            RankState& rs = global.ranks[r];
            const size_t n = rs.nodes.size();
            if (n < 4) continue;
            // Proxy wall time per node: MPI nodes use their own midpoint;
            // COMPUTE gaps inherit the last known wall time (gap start).
            std::vector<double> wt(n);
            double last_w = -DBL_MAX;
            for (size_t i = 0; i < n; i++) {
                if (rs.nodes[i].wall_time_sec >= 0.0)
                    last_w = rs.nodes[i].wall_time_sec;
                wt[i] = last_w;
            }
            if (last_w <= -DBL_MAX) continue;
            // Anchor: the LAMMPS Timer sync mode brackets the run loop
            // with MPI_Barrier (barrier_start / barrier_stop) — exactly
            // the endpoints of "Loop time" in log.lammps. The single
            // captured step therefore spans [last BARRIER before the last
            // P2P node, first BARRIER after it]; the setup prefix and the
            // post-step thermo printout (ALLREDUCE burst + BCAST + COMPUTE)
            // both lie outside this pair and are trimmed, while the step's
            // halo P2P stays complete.
            size_t p_last = n;   // last P2P node (thermo printout has none)
            for (size_t i = n; i-- > 0; )
                if (is_p2p_type(rs.nodes[i].type)) { p_last = i; break; }
            size_t b_tail = n, b_head = n;
            if (p_last < n) {
                for (size_t i = p_last + 1; i < n; i++)
                    if (rs.nodes[i].type == COMM_BARRIER) { b_tail = i; break; }
                for (size_t i = p_last; i-- > 0; )
                    if (rs.nodes[i].type == COMM_BARRIER) { b_head = i; break; }
            }
            double end_w = -DBL_MAX;
            double start_w = -DBL_MAX;
            if (b_head < n && b_tail < n) {
                double dur = rs.nodes[b_tail].wall_duration_sec > 0.0
                             ? rs.nodes[b_tail].wall_duration_sec : 0.0;
                start_w = rs.nodes[b_head].wall_time_sec;
                end_w = rs.nodes[b_tail].wall_time_sec + dur / 2.0;
                fprintf(stderr,
                        "Trim setup: rank %d step anchored on BARRIER pair "
                        "(span %.1f us vs Loop time %.1f us)\n",
                        r, (end_w - start_w) * 1e6, loop_time * 1e6);
            } else {
                // Fallback: end of loop = wall stop of the last in-loop
                // call (the post-run printout suffix lies after the Loop
                // timer stops and is skipped).
                size_t anchor = n;
                for (size_t i = n; i-- > 0; ) {
                    NodeType t = rs.nodes[i].type;
                    if (t == COMPUTE || t == COMM_BCAST || t == COMM_BARRIER ||
                        t == COMM_REDUCE || t == COMM_GATHER) continue;
                    if (rs.nodes[i].wall_time_sec < 0.0) continue;
                    anchor = i;
                    break;
                }
                if (anchor < n) {
                    double dur = rs.nodes[anchor].wall_duration_sec > 0.0
                                 ? rs.nodes[anchor].wall_duration_sec : 0.0;
                    end_w = rs.nodes[anchor].wall_time_sec + dur / 2.0;
                } else {
                    for (size_t i = 0; i < n; i++) {
                        if (rs.nodes[i].wall_time_sec >= 0.0) {
                            double dur = rs.nodes[i].wall_duration_sec > 0.0
                                         ? rs.nodes[i].wall_duration_sec : 0.0;
                            end_w = std::max(end_w, rs.nodes[i].wall_time_sec + dur / 2.0);
                        }
                    }
                }
                if (end_w <= -DBL_MAX) continue;
                start_w = end_w - loop_time;
            }
            if (end_w <= -DBL_MAX || start_w <= -DBL_MAX) continue;
            // Fallback rollback: keep the last P2P burst complete when
            // the fallback Loop-time anchor would cut into the step's
            // halo exchange (trace wall timestamps diverging from the
            // LAMMPS Loop time). No-op for the BARRIER-anchored path
            // since the step-start BARRIER precedes the burst.
            if (p_last < n) {
                size_t burst_start = p_last;
                while (burst_start > 0) {
                    NodeType t = rs.nodes[burst_start - 1].type;
                    if (is_p2p_type(t) || t == COMPUTE) burst_start--;
                    else break;
                }
                double burst_w = wt[burst_start];
                if (burst_w >= 0.0 && burst_w < start_w) {
                    fprintf(stderr,
                            "Trim setup: rank %d window start rolled back %.1f us "
                            "(node %zu) to keep last P2P burst\n",
                            r, (start_w - burst_w) * 1e6, burst_start);
                    start_w = burst_w;
                }
            }
            size_t i0 = 0, i1 = n; // keep [i0, i1)
            while (i0 < n && wt[i0] < start_w) i0++;
            while (i1 > i0 && wt[i1 - 1] > end_w) i1--;
            if (i0 == 0 && i1 == n) continue;
            if (i1 - i0 < 4) {
                fprintf(stderr, "Trim setup: rank %d window degenerate, skipping\n", r);
                continue;
            }
            // Clip the boundary COMPUTE gap by its wall overlap with the
            // window instead of zeroing it (coarse-call traces like few-rank
            // runs would otherwise lose a whole in-loop compute gap).
            if (rs.nodes[i0].type == COMPUTE) {
                double gstart = wt[i0];
                double gend = end_w;
                for (size_t j = i0 + 1; j < i1; j++) {
                    if (rs.nodes[j].wall_time_sec >= 0.0) {
                        gend = rs.nodes[j].wall_time_sec;
                        break;
                    }
                }
                if (gend > gstart) {
                    double lo = std::max(gstart, start_w);
                    double hi = std::min(gend, end_w);
                    double frac = (hi > lo) ? (hi - lo) / (gend - gstart) : 0.0;
                    if (frac > 1.0) frac = 1.0;
                    rs.nodes[i0].compute_time_sec *= frac;
                    rs.nodes[i0].compute_cycles *= frac;
                    rs.nodes[i0].compute_ops *= frac;
                }
            }
            size_t removed = n - (i1 - i0);
            rs.nodes.erase(rs.nodes.begin() + i1, rs.nodes.end());
            rs.nodes.erase(rs.nodes.begin(), rs.nodes.begin() + i0);
            rs.nodes.front().predecessors.clear();
            window_total += removed;
            fprintf(stderr, "Trim setup: rank %d window removed %zu nodes (head %zu / tail %zu, %zu remain)\n",
                    r, removed, i0, n - i1, rs.nodes.size());
            // Report the retained compute volume against the Loop time so
            // that residual setup compute or over-trimming is visible
            // (SimGrid frequency calibration will flag outliers).
            double win_cycles = 0.0;
            for (auto& nd : rs.nodes)
                if (nd.type == COMPUTE) win_cycles += nd.compute_cycles;
            fprintf(stderr,
                    "Trim setup: rank %d window compute = %.1f us @%.1fGHz (Loop time %.1f us, ratio %.2f)\n",
                    r, win_cycles / g_cpu_freq_ghz / 1e9 * 1e6, g_cpu_freq_ghz, loop_time * 1e6,
                    win_cycles / g_cpu_freq_ghz / 1e9 / loop_time);
        }
        if (window_total > 0) {
            renumber_ids(global);
            fprintf(stderr, "Trim setup: loop-time window removed %zu nodes in total\n", window_total);
        }
    } else {
        fprintf(stderr, "Trim setup: no LAMMPS log Loop time found, skipping window trim\n");
    }

    if (trimmed_total > 0) {
        renumber_ids(global);
        fprintf(stderr, "Trim setup: removed %zu nodes in total\n", trimmed_total);
    }
}

/* Burst compaction (env CCDG_COMPACT): fold every maximal run of
 * {COMPUTE, P2P} nodes into aggregate nodes:
 *   COMPUTE          -> 1 aggregate COMPUTE (cycles/time summed)
 *   SEND/ISEND       -> 1 aggregate SEND per comm_dst (bytes summed)
 *   RECV/IRECV       -> 1 aggregate IRECV per comm_src (bytes summed)
 *   WAIT/WAITANY/ALL -> 1 aggregate WAITALL
 * Cross edges are remapped to the aggregate nodes and de-duplicated;
 * total comm bytes must be exactly conserved. Run after
 * build_cross_edges + eliminate_cycles so the validated SEND<->WAIT
 * matching is reused. */
static void compact_bursts(GlobalState& global) {
    uint64_t bytes_before = 0;
    size_t nodes_before = 0;
    for (auto& rs : global.ranks) {
        nodes_before += rs.nodes.size();
        for (auto& n : rs.nodes) bytes_before += n.comm_bytes;
    }

    std::map<uint64_t, uint64_t> idmap;  // old node id -> aggregate node id
    size_t folded_bursts = 0;

    for (int r = 0; r < global.num_ranks; r++) {
        RankState& rs = global.ranks[r];
        std::vector<CCDGNode> out;
        size_t i = 0;
        while (i < rs.nodes.size()) {
            if (!is_burst_type(rs.nodes[i].type)) {
                idmap[rs.nodes[i].id] = rs.nodes[i].id;
                out.push_back(rs.nodes[i]);
                i++;
                continue;
            }
            size_t j = i;
            while (j < rs.nodes.size() && is_burst_type(rs.nodes[j].type)) j++;

            // ---- fold burst [i, j) ----
            CCDGNode ac; bool has_ac = false;
            std::map<int, CCDGNode> sends;  // comm_dst -> aggregate SEND
            std::map<int, CCDGNode> recvs;  // comm_src -> aggregate IRECV
            CCDGNode aw; bool has_aw = false;
            for (size_t t = i; t < j; t++) {
                const CCDGNode& n = rs.nodes[t];
                switch (n.type) {
                case COMPUTE:
                    if (!has_ac) { ac = n; has_ac = true; }
                    else {
                        ac.compute_cycles += n.compute_cycles;
                        ac.compute_time_sec += n.compute_time_sec;
                        ac.compute_ops += n.compute_ops;
                    }
                    break;
                case COMM_SEND: case COMM_ISEND: {
                    auto it = sends.find(n.comm_dst);
                    if (it == sends.end()) {
                        CCDGNode s = n; s.type = COMM_SEND; s.pending_req_id = -1;
                        sends[n.comm_dst] = s;
                    } else {
                        it->second.comm_bytes += n.comm_bytes;
                        it->second.comm_count += n.comm_count;
                    }
                    break;
                }
                case COMM_RECV: case COMM_IRECV: {
                    auto it = recvs.find(n.comm_src);
                    if (it == recvs.end()) {
                        CCDGNode v = n; v.type = COMM_IRECV; v.pending_req_id = -1;
                        recvs[n.comm_src] = v;
                    } else {
                        it->second.comm_bytes += n.comm_bytes;
                        it->second.comm_count += n.comm_count;
                    }
                    break;
                }
                default: {  // WAIT / WAITANY / WAITALL
                    double wd = n.wall_duration_sec > 0 ? n.wall_duration_sec : 0;
                    if (!has_aw) {
                        aw = n; aw.type = COMM_WAITALL; aw.pending_req_id = -1;
                        aw.wall_duration_sec = wd; has_aw = true;
                    } else aw.wall_duration_sec += wd;
                    break;
                }
                }
            }

            std::vector<CCDGNode> repl;
            if (has_ac) { ac.id = global.global_node_counter++; ac.predecessors.clear(); repl.push_back(ac); }
            for (auto& kv : recvs) { CCDGNode v = kv.second; v.id = global.global_node_counter++; v.predecessors.clear(); repl.push_back(v); }
            for (auto& kv : sends) { CCDGNode s = kv.second; s.id = global.global_node_counter++; s.predecessors.clear(); repl.push_back(s); }
            if (has_aw) { aw.id = global.global_node_counter++; aw.predecessors.clear(); repl.push_back(aw); }
            if (repl.empty()) { i = j; continue; }

            // old id -> aggregate id mapping (same map iteration order as repl)
            size_t pos = 0;
            uint64_t id_ac = has_ac ? repl[pos++].id : 0;
            std::map<int, uint64_t> rid_recv, rid_send;
            for (auto& kv : recvs) rid_recv[kv.first] = repl[pos++].id;
            for (auto& kv : sends) rid_send[kv.first] = repl[pos++].id;
            uint64_t id_aw = has_aw ? repl[pos++].id : 0;
            for (size_t t = i; t < j; t++) {
                const CCDGNode& n = rs.nodes[t];
                uint64_t target;
                switch (n.type) {
                case COMPUTE: target = has_ac ? id_ac : repl.front().id; break;
                case COMM_SEND: case COMM_ISEND: target = rid_send[n.comm_dst]; break;
                case COMM_RECV: case COMM_IRECV: target = rid_recv[n.comm_src]; break;
                default: target = has_aw ? id_aw : repl.back().id; break;
                }
                idmap[n.id] = target;
            }

            if (j - i > repl.size()) folded_bursts++;
            for (auto& rn : repl) out.push_back(rn);
            i = j;
        }
        // rebuild the intra-rank program-order chain
        for (size_t k = 0; k < out.size(); k++) {
            out[k].predecessors.clear();
            if (k > 0) out[k].predecessors.push_back(out[k - 1].id);
        }
        rs.nodes.swap(out);
    }

    // remap + dedup cross edges
    std::vector<std::pair<uint64_t, uint64_t>> new_edges;
    std::set<std::pair<uint64_t, uint64_t>> seen;
    size_t dropped = 0;
    for (auto& e : global.cross_edges) {
        uint64_t u = idmap[e.first], v = idmap[e.second];
        if (u == v || !seen.insert({u, v}).second) { dropped++; continue; }
        new_edges.push_back({u, v});
    }
    global.cross_edges.swap(new_edges);

    renumber_ids(global);

    uint64_t bytes_after = 0;
    size_t nodes_after = 0;
    for (auto& rs : global.ranks) {
        nodes_after += rs.nodes.size();
        for (auto& n : rs.nodes) bytes_after += n.comm_bytes;
    }
    fprintf(stderr, "Compact bursts: nodes %zu -> %zu (%.1fx), bursts folded %zu, cross edges dedup -%zu\n",
            nodes_before, nodes_after,
            (double)nodes_before / (double)std::max<size_t>(nodes_after, 1),
            folded_bursts, dropped);
    if (bytes_after != bytes_before)
        fprintf(stderr, "WARNING: comm bytes NOT conserved: %lu -> %lu\n",
                (unsigned long)bytes_before, (unsigned long)bytes_after);
    else
        fprintf(stderr, "Compact bursts: comm bytes conserved: %lu\n",
                (unsigned long)bytes_before);
}

static int on_allreduce_cb(const dumpi_allreduce* prm, uint16_t thread,
                           const dumpi_time* cpu, const dumpi_time* wall,
                           const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->count * compute_datatype_size(prm->datatype);
    rs.add_collective_node(COMM_ALLREDUCE, s, bytes * ctx->global->num_ranks);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_barrier_cb(const dumpi_barrier* prm, uint16_t thread,
                         const dumpi_time* cpu, const dumpi_time* wall,
                         const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    rs.add_collective_node(COMM_BARRIER, s);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_bcast_cb(const dumpi_bcast* prm, uint16_t thread,
                       const dumpi_time* cpu, const dumpi_time* wall,
                       const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->count * compute_datatype_size(prm->datatype);
    rs.add_collective_node(COMM_BCAST, s, bytes * ctx->global->num_ranks, prm->root);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_gather_cb(const dumpi_gather* prm, uint16_t thread,
                        const dumpi_time* cpu, const dumpi_time* wall,
                        const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->sendcount * compute_datatype_size(prm->sendtype);
    rs.add_collective_node(COMM_GATHER, s, bytes, prm->root);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_allgather_cb(const dumpi_allgather* prm, uint16_t thread,
                           const dumpi_time* cpu, const dumpi_time* wall,
                           const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->sendcount * compute_datatype_size(prm->sendtype);
    rs.add_collective_node(COMM_ALLGATHER, s, bytes * ctx->global->num_ranks);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_scatter_cb(const dumpi_scatter* prm, uint16_t thread,
                         const dumpi_time* cpu, const dumpi_time* wall,
                         const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->sendcount * compute_datatype_size(prm->sendtype);
    rs.add_collective_node(COMM_SCATTER, s, bytes, prm->root);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_alltoall_cb(const dumpi_alltoall* prm, uint16_t thread,
                          const dumpi_time* cpu, const dumpi_time* wall,
                          const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->sendcount * compute_datatype_size(prm->sendtype);
    rs.add_collective_node(COMM_ALLTOALL, s, bytes * ctx->global->num_ranks);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

static int on_reduce_cb(const dumpi_reduce* prm, uint16_t thread,
                        const dumpi_time* cpu, const dumpi_time* wall,
                        const dumpi_perfinfo* perf, void* userarg) {
    auto* ctx = static_cast<CallbackContext*>(userarg);
    RankState& rs = ctx->global->ranks[ctx->rank];
    double s, e;
    CPU_GAP_BEGIN(rs, cpu, s, e);
    uint64_t bytes = prm->count * compute_datatype_size(prm->datatype);
    rs.add_collective_node(COMM_REDUCE, s, bytes, prm->root);
    rs.set_wall_time(wall_time_to_sec(wall), wall_duration(wall));
    CPU_GAP_END(rs, e);
    return 0;
}

/* ================================================================
 * Parse a single rank's trace file
 * ================================================================ */
static bool parse_rank_trace(GlobalState& global, int rank,
                             const std::string& binfile,
                             const std::string& metafile) {
    CallbackContext ctx;
    ctx.global = &global;
    ctx.rank = rank;
    // Set the global counter pointer for this rank
    global.ranks[rank].global_counter = &global.global_node_counter;

    dumpi_profile* profile = undumpi_open(binfile.c_str());
    if (!profile) {
        fprintf(stderr, "ERROR: undumpi_open failed for %s\n", binfile.c_str());
        return false;
    }

    // Read header
    dumpi_header* header = undumpi_read_header(profile);
    if (header) {
        dumpi_free_header(header);
    }

    // Set up callbacks
    libundumpi_callbacks cb;
    libundumpi_clear_callbacks(&cb);
    cb.on_send      = on_send_cb;
    cb.on_recv      = on_recv_cb;
    cb.on_isend     = on_isend_cb;
    cb.on_irecv     = on_irecv_cb;
    cb.on_wait      = on_wait_cb;
    cb.on_waitany   = on_waitany_cb;
    cb.on_waitall   = on_waitall_cb;
    cb.on_sendrecv  = on_sendrecv_cb;
    cb.on_allreduce = on_allreduce_cb;
    cb.on_barrier   = on_barrier_cb;
    cb.on_bcast     = on_bcast_cb;
    cb.on_gather    = on_gather_cb;
    cb.on_allgather = on_allgather_cb;
    cb.on_scatter   = on_scatter_cb;
    cb.on_alltoall  = on_alltoall_cb;
    cb.on_reduce    = on_reduce_cb;

    // Read stream
    int ret = undumpi_read_stream(profile, &cb, &ctx, false);
    if (!ret) {
        fprintf(stderr, "WARNING: undumpi_read_stream returned %d for %s\n",
                ret, binfile.c_str());
    }

    // Read footer
    dumpi_footer* footer = undumpi_read_footer(profile);
    if (footer) {
        dumpi_free_footer(footer);
    }

    undumpi_close(profile);
    return true;
}

/* ================================================================
 * Find trace files for all ranks
 * ================================================================ */
static bool find_trace_files(const std::string& trace_dir,
                             std::vector<std::string>& bin_files,
                             std::string& meta_file,
                             int& num_ranks) {
    // Find meta file
    glob_t g;
    std::string meta_pattern = trace_dir + "/dumpi-*.meta";
    int ret = glob(meta_pattern.c_str(), 0, nullptr, &g);
    if (ret != 0 || g.gl_pathc == 0) {
        fprintf(stderr, "ERROR: no .meta files found in %s\n", trace_dir.c_str());
        globfree(&g);
        return false;
    }
    meta_file = g.gl_pathv[0];
    globfree(&g);

    // Find bin files
    std::string bin_pattern = trace_dir + "/dumpi-*.bin";
    ret = glob(bin_pattern.c_str(), 0, nullptr, &g);
    if (ret != 0 || g.gl_pathc == 0) {
        fprintf(stderr, "ERROR: no .bin files found in %s\n", trace_dir.c_str());
        globfree(&g);
        return false;
    }

    // Sort by rank number (files are named dumpi-<timestamp>-<rank>.bin)
    for (size_t i = 0; i < g.gl_pathc; i++) {
        bin_files.push_back(g.gl_pathv[i]);
    }
    std::sort(bin_files.begin(), bin_files.end());
    num_ranks = (int)bin_files.size();
    globfree(&g);
    return true;
}

/* ================================================================
 * JSON output
 * ================================================================ */
static void print_json(const GlobalState& global) {
    printf("{\n");
    printf("  \"num_ranks\": %d,\n", global.num_ranks);
    printf("  \"nodes\": [\n");

    bool first_node = true;
    for (int r = 0; r < global.num_ranks; r++) {
        const RankState& rs = global.ranks[r];
        for (const auto& node : rs.nodes) {
            if (!first_node) printf(",\n");
            first_node = false;
            printf("    {\n");
            printf("      \"id\": %lu,\n", node.id);
            printf("      \"rank\": %d,\n", node.rank);
            printf("      \"type\": \"%s\"", node_type_str(node.type));
            if (node.compute_cycles > 0)
                printf(",\n      \"compute_cycles\": %.0f", node.compute_cycles);
            if (node.compute_time_sec > 0)
                printf(",\n      \"compute_time_sec\": %.9f", node.compute_time_sec);
            if (node.compute_ops > 0)
                printf(",\n      \"compute_ops\": %.0f", node.compute_ops);
            if (node.wall_time_sec >= 0)
                printf(",\n      \"wall_time_sec\": %.9f", node.wall_time_sec);
            if (node.wall_duration_sec >= 0)
                printf(",\n      \"wall_duration_sec\": %.9f", node.wall_duration_sec);
            if (node.type != COMPUTE) {
                if (node.comm_src >= 0)
                    printf(",\n      \"comm_src\": %d", node.comm_src);
                if (node.comm_dst >= 0)
                    printf(",\n      \"comm_dst\": %d", node.comm_dst);
                if (node.comm_tag >= 0)
                    printf(",\n      \"comm_tag\": %d", node.comm_tag);
                if (node.comm_bytes > 0)
                    printf(",\n      \"comm_bytes\": %lu", node.comm_bytes);
                if (node.comm_count > 0)
                    printf(",\n      \"comm_count\": %d", node.comm_count);
                if (node.collective_root >= 0)
                    printf(",\n      \"collective_root\": %d", node.collective_root);
                if (node.pending_req_id >= 0)
                    printf(",\n      \"pending_req_id\": %d", node.pending_req_id);
            }
            // Predecessors (intra-rank)
            if (!node.predecessors.empty()) {
                printf(",\n      \"predecessors\": [");
                for (size_t i = 0; i < node.predecessors.size(); i++) {
                    if (i > 0) printf(", ");
                    printf("%lu", node.predecessors[i]);
                }
                printf("]");
            }
            printf("\n    }");
        }
    }
    printf("\n  ],\n");

    // Cross-rank edges
    printf("  \"cross_rank_edges\": [\n");
    bool first_edge = true;
    for (const auto& edge : global.cross_edges) {
        if (!first_edge) printf(",\n");
        first_edge = false;
        printf("    { \"src_node\": %lu, \"dst_node\": %lu }", edge.first, edge.second);
    }
    printf("\n  ]\n");
    printf("}\n");
}

/* ================================================================
 * Main
 * ================================================================ */
int main(int argc, char** argv) {
    // CPU-frequency assumption and legacy midpoint mode are read from the
    // environment before any trace parsing (CCDG_CPU_FREQ_GHZ default 2.5;
    // CCDG_LEGACY_CPU_MID=1 reproduces the old (start+stop)/2 extraction).
    if (getenv("CCDG_CPU_FREQ_GHZ")) {
        g_cpu_freq_ghz = atof(getenv("CCDG_CPU_FREQ_GHZ"));
        fprintf(stderr, "CCDG_CPU_FREQ_GHZ = %.3f\n", g_cpu_freq_ghz);
    }
    if (getenv("CCDG_LEGACY_CPU_MID")) {
        g_legacy_cpu_mid = true;
        fprintf(stderr, "Legacy CPU midpoint mode enabled (CCDG_LEGACY_CPU_MID)\n");
    }
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <trace_dir>\n", argv[0]);
        fprintf(stderr, "  trace_dir: directory containing dumpi-*.bin and dumpi-*.meta files\n");
        return 1;
    }

    std::string trace_dir = argv[1];
    std::vector<std::string> bin_files;
    std::string meta_file;
    int num_ranks = 0;

    if (!find_trace_files(trace_dir, bin_files, meta_file, num_ranks)) {
        return 1;
    }

    fprintf(stderr, "Found %d rank trace files in %s\n", num_ranks, trace_dir.c_str());

    GlobalState global(num_ranks);

    // Parse each rank's trace
    for (int r = 0; r < num_ranks; r++) {
        fprintf(stderr, "Parsing rank %d/%d: %s\n", r, num_ranks-1, bin_files[r].c_str());
        if (!parse_rank_trace(global, r, bin_files[r], meta_file)) {
            fprintf(stderr, "ERROR: failed to parse rank %d\n", r);
            return 1;
        }
        fprintf(stderr, "  -> %zu nodes\n", global.ranks[r].nodes.size());
    }

    // Print statistics
    size_t total_nodes = 0;
    size_t total_compute = 0;
    size_t total_comm = 0;
    double total_compute_time = 0;
    for (int r = 0; r < num_ranks; r++) {
        const RankState& rs = global.ranks[r];
        total_nodes += rs.nodes.size();
        for (const auto& n : rs.nodes) {
            if (n.type == COMPUTE) {
                total_compute++;
                total_compute_time += n.compute_time_sec;
            } else {
                total_comm++;
            }
        }
    }
    // Optional: trim the setup prefix before matching (CCDG_TRIM_SETUP)
    if (getenv("CCDG_TRIM_SETUP")) {
        fprintf(stderr, "Setup trimming enabled (CCDG_TRIM_SETUP)\n");
        trim_setup(global, trace_dir);
    }

    // Build cross-rank edges from matched send/recv-wait pairs
    build_cross_edges(global);

    // Drop cross edges that would deadlock under collective-sync semantics
    eliminate_cycles(global);

    // Optional: fold P2P bursts into aggregate nodes (CCDG_COMPACT)
    if (getenv("CCDG_COMPACT")) {
        fprintf(stderr, "Burst compaction enabled (CCDG_COMPACT)\n");
        compact_bursts(global);
    }

    // Recompute statistics after optional trimming/compaction
    total_nodes = 0;
    total_compute = 0;
    total_comm = 0;
    total_compute_time = 0;
    for (int r = 0; r < num_ranks; r++) {
        const RankState& rs = global.ranks[r];
        total_nodes += rs.nodes.size();
        for (const auto& n : rs.nodes) {
            if (n.type == COMPUTE) {
                total_compute++;
                total_compute_time += n.compute_time_sec;
            } else {
                total_comm++;
            }
        }
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Statistics:\n");
    fprintf(stderr, "  Total nodes: %zu\n", total_nodes);
    fprintf(stderr, "  Compute nodes: %zu\n", total_compute);
    fprintf(stderr, "  Communication nodes: %zu\n", total_comm);
    fprintf(stderr, "  Total compute time: %.3f seconds\n", total_compute_time);
    fprintf(stderr, "  Cross-rank edges: %zu\n", global.cross_edges.size());
    fprintf(stderr, "\n");

    // Output JSON
    print_json(global);

    return 0;
}