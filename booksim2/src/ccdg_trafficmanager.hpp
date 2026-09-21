// $Id$

/*
 Copyright (c) 2007-2015, Trustees of The Leland Stanford Junior University
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice, this
 list of conditions and the following disclaimer in the documentation and/or
 other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _CCDG_TRAFFICMANAGER_HPP_
#define _CCDG_TRAFFICMANAGER_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstdint>

#include "config_utils.hpp"
#include "stats.hpp"
#include "trafficmanager.hpp"
#include "flit.hpp"

class CCDGTrafficManager : public TrafficManager {

protected:

  // ============ CCDG Data Structures ============

  struct CCDGNodeInfo {
    uint64_t    id;
    int         rank;
    std::string type;
    double      compute_cycles;
    double      compute_ops;    // ops estimate from the trace (v3; 0 for legacy CCDG)
    int         comm_src;
    int         comm_dst;
    int         comm_tag;
    uint64_t    comm_bytes;
    int         coll_root;      // root rank for rooted collectives (-1 if none)
    int         wse_wavefront_idx; // manager-level replay metadata (-1 for CCDG)
    int         wse_stage_idx;     // flattened compiler stage index
    std::string wse_kind;          // "command" or "branch"
    std::string wse_mode;          // "multicast" or "reduction"
    // Compile-time earliest start time (NoC cycles) from the deterministic
    // schedule (ccdg_schedule_file); 0 = unconstrained (responsive)
    double      sched_est;
    // Indices into _cross_edges where this node is the destination
    std::vector<size_t> dep_edge_indices;
    // Indices into _cross_edges where this node is the source
    std::vector<size_t> src_edge_indices;
  };

  struct CrossRankEdgeInfo {
    uint64_t src_node_id;
    uint64_t dst_node_id;
    int      src_rank;
    int      dst_rank;
    int      msg_id;       // assigned during simulation, -1 if not yet assigned
    bool     resolved;
  };

  // ============ PE State ============

  enum PEStateEnum { PE_COMPUTE, PE_BLOCKED, PE_GATED, PE_BACKPRESSURE, PE_DONE };

  struct PERankState {
    size_t     current_node_idx;   // index into _rank_nodes[rank]
    double     remaining_cycles;
    PEStateEnum state;
    // Set of cross-edge indices that this rank is waiting on
    std::set<size_t> pending_edges;
    int        coll_seq;           // number of collective nodes entered so far
    // Statistics (accumulated over one _SingleSim run)
    double     blocked_cycles;     // cycles spent in PE_BLOCKED (exposed comm wait)
    double     congestion_cycles;     // cycles spent in PE_BACKPRESSURE (exposed comm wait)
    double     compute_cycles_acc; // cycles spent in PE_COMPUTE (roofline dwell)
    double     sched_wait_cycles;  // orchestration slots: parked before the
                                   // compile-time EST release time
  };

  // Per-rank state for collective traffic accounting (non-blocking model:
  // collectives inject their expanded packets and advance immediately;
  // the network cost shows up through contention and completion drain)
  struct PECollState {
    int pending_out;     // issued collective packets not yet retired
    int injected;        // expanded packets already injected at the current
                         // collective node (partial progress under backpressure)
    PECollState() : pending_out(0), injected(0) {}
  };

  // Bookkeeping for an in-flight collective packet
  struct CollMsgInfo {
    int src;
    int dst;
  };

  // ============ CCDG Data ============

  std::vector<std::vector<CCDGNodeInfo> > _rank_nodes;
  std::vector<CrossRankEdgeInfo> _cross_edges;
  std::vector<PERankState> _pe_state;

  // ============ Configuration ============

  std::string _ccdg_file;
  double      _cpu_freq_ghz;      // CPU frequency in GHz
  double      _noc_freq_ghz;      // NoC frequency in GHz
  double      _freq_ratio;        // CPU_freq / NoC_freq
  double      _ccdg_compute_rate; // PE compute rate in ops/cycle (roofline);
                                  // default = _freq_ratio (legacy-compatible)
  double      _ccdg_compute_capability; // PE compute capability in ops/s (input param;
                                        // 0 = unset, fall back to compute_rate)
  int         _ccdg_inject_queue_depth; // injection queue depth in flits for PE
                                        // backpressure (0 = auto: max single-packet
                                        // flit count; queue-full stalls the PE)
  uint64_t    _bp_check_hits;   // debug: backpressure check triggered count
  uint64_t    _bp_state_cycles; // debug: cycles spent in PE_BACKPRESSURE
  int         _bp_max_seen;     // debug: max queue size observed at trigger
  int         _flit_size_bytes;
  int         _num_ranks;

  // Deterministic orchestration (prescriptive schedule): optional EST table
  // from ccdg_schedule_file ("node_id est_cycles" per line). When enabled the
  // PE parks at each node until its release time (sched_wait), executing the
  // compile-time schedule instead of responsive dependency resolution.
  bool        _sched_enabled;
  std::string _ccdg_schedule_file;

  // ============ Message Tracking ============

  int  _msg_id_counter;
  // Map from msg_id to list of cross-edge indices
  std::map<int, std::vector<size_t> > _msg_id_to_edge_idx;

  // ============ Collective Tracking ============

  std::vector<PECollState> _coll_state;
  // msg_id -> (src, dst, seq) for collective packets
  std::map<int, CollMsgInfo> _coll_msg_info;
  // barrier seq -> set of ranks that have arrived
  std::map<int, std::set<int> > _barrier_arrived;

  // ============ Statistics ============

  int _total_sim_cycles;
  int _total_packets_sent;
  int _total_packets_received;
  uint64_t _injected_flits_total;
  uint64_t _injected_packets_total;
  uint64_t _injection_capacity_slots;
  uint64_t _backlogged_injection_slots;
  uint64_t _flit_queue_cycles_total;
  uint64_t _packet_queue_cycles_total;

  // ============ Internal Methods ============

  // Parse a CCDG JSON file
  bool _parseCCDG(const std::string &filename);

  // Advance a PE to the next node(s) in its CCDG sequence
  void _advancePE(int rank);

  // Generate a packet for a communication node
  void _issueCCDGPacket(int rank, const CCDGNodeInfo &node);

  // Inject one packet of flits from src to dst carrying msg_id
  void _injectPacket(int src, int dst, uint64_t bytes, int msg_id);

  // Hooks used by WSETrafficManager without duplicating the network engine.
  virtual void _OnPacketIssued(int msg_id, const CCDGNodeInfo &node);
  virtual void _OnPacketRetired(int msg_id, int dest);

  // Compute the expanded send destinations and expected receive count
  // for a collective operation on a given rank
  void _collectivePattern(const std::string &type, int rank, int root,
                          std::vector<int> &dests, int &expected_in) const;

  // Get the number of flits for a given message size
  int _getFlitCount(uint64_t bytes) const;

  // CCDG-specific aggregate metrics.  A slot is one node/subnetwork injection
  // opportunity in one cycle; a backlogged slot has a matching queued flit.
  double _AverageFlitQueueCycles() const;
  double _AveragePacketQueueCycles() const;
  double _AverageInjectionRate() const;
  double _InjectionSaturationRatio() const;
  double _SaturatedInjectionRate() const;
  double _ExposedCommunicationToComputeRatio() const;
  double _CommunicationToComputeRatio() const;

  // Load the EST schedule file and attach each node's release time
  bool _loadScheduleFile();

  // ============ Overridden Virtual Methods ============

  virtual void _RetireFlit(Flit *f, int dest);
  virtual bool _SingleSim();
  virtual void _ClearStats();
  virtual void _UpdateOverallStats();
  virtual std::string _OverallStatsCSV(int c = 0) const;

public:

  CCDGTrafficManager(const Configuration &config, const std::vector<Network *> &net);
  virtual ~CCDGTrafficManager();

  virtual void WriteStats(std::ostream &os = std::cout) const;
  virtual void DisplayStats(std::ostream &os = std::cout) const;
  virtual void DisplayOverallStats(std::ostream &os = std::cout) const;

};

#endif
