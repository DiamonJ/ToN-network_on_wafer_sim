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

#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>

#include "booksim.hpp"
#include "booksim_config.hpp"
#include "ccdg_trafficmanager.hpp"
#include "random_utils.hpp"
#include "vc.hpp"
#include "packet_reply_info.hpp"

CCDGTrafficManager::CCDGTrafficManager(const Configuration &config,
                                       const vector<Network *> &net)
  : TrafficManager(config, net),
    _msg_id_counter(0), _total_sim_cycles(0),
    _total_packets_sent(0), _total_packets_received(0)
{
  // Read CCDG-specific configuration
  _ccdg_file = config.GetStr("ccdg_file");
  _cpu_freq_ghz = config.GetFloat("cpu_frequency_ghz");
  if (_cpu_freq_ghz <= 0.0) _cpu_freq_ghz = 2.5; // default 2.5 GHz
  _noc_freq_ghz = config.GetFloat("noc_frequency_ghz");
  if (_noc_freq_ghz <= 0.0) _noc_freq_ghz = 1.0; // default 1.0 GHz
  _freq_ratio = _cpu_freq_ghz / _noc_freq_ghz;
  // Input parameter: compute capability in ops/s (target hardware).
  // compute time (sec) = compute_ops / capability; cycles = time x noc_freq.
  // Equivalent per-cycle rate = capability / (noc_freq_ghz * 1e9) ops/cycle.
  // 0 = unset: fall back to ccdg_compute_rate (ops/cycle), then _freq_ratio (legacy).
  _ccdg_compute_capability = config.GetFloat("ccdg_compute_capability");
  _ccdg_compute_rate = config.GetFloat("ccdg_compute_rate");
  if (_ccdg_compute_capability > 0.0) {
    _ccdg_compute_rate = _ccdg_compute_capability / (_noc_freq_ghz * 1e9);
  }
  if (_ccdg_compute_rate <= 0.0) _ccdg_compute_rate = _freq_ratio; // legacy default
  _ccdg_inject_queue_depth = config.GetInt("ccdg_inject_queue_depth");
  _bp_check_hits = 0;
  _bp_state_cycles = 0;
  _bp_max_seen = 0;
  _flit_size_bytes = config.GetInt("flit_size_bytes");
  if (_flit_size_bytes <= 0) _flit_size_bytes = 8; // default 8 bytes per flit
  _ccdg_schedule_file = config.GetStr("ccdg_schedule_file");
  _sched_enabled = !_ccdg_schedule_file.empty() && _ccdg_schedule_file != "none";

  cout << "CCDGTrafficManager: loading CCDG from " << _ccdg_file << endl;
  cout << "  CPU frequency: " << _cpu_freq_ghz << " GHz" << endl;
  cout << "  NoC frequency: " << _noc_freq_ghz << " GHz" << endl;
  cout << "  Frequency ratio: " << _freq_ratio << endl;
  if (_ccdg_compute_capability > 0.0) {
    cout << "  Compute capability: " << _ccdg_compute_capability
         << " ops/s (" << _ccdg_compute_rate << " ops/cycle @ "
         << _noc_freq_ghz << " GHz)" << endl;
  } else {
    cout << "  Compute rate: " << _ccdg_compute_rate << " ops/cycle (roofline)" << endl;
  }
  cout << "  Flit size: " << _flit_size_bytes << " bytes" << endl;
  if (_sched_enabled) {
    cout << "  Deterministic orchestration: schedule file "
         << _ccdg_schedule_file << endl;
  } else {
    cout << "  Deterministic orchestration: disabled (responsive execution)" << endl;
  }

  // Parse the CCDG file
  if (!_parseCCDG(_ccdg_file)) {
    Error("Failed to parse CCDG file: " + _ccdg_file);
  }

  if (_sched_enabled && !_loadScheduleFile()) {
    Error("Failed to load schedule file: " + _ccdg_schedule_file);
  }

  _num_ranks = (int)_rank_nodes.size();
  cout << "  Loaded " << _num_ranks << " ranks, "
       << _cross_edges.size() << " cross-rank edges" << endl;

  // Injection queue depth for PE backpressure: 0 = auto = max single-packet
  // flit count. When the queue cannot hold the next whole packet the PE
  // stalls (PE_BACKPRESSURE), so network congestion propagates back into
  // the CCDG schedule (real-hardware credit backpressure).
  if (_ccdg_inject_queue_depth <= 0) {
    uint64_t max_bytes = 0;
    for (int r = 0; r < _num_ranks; r++) {
      for (size_t i = 0; i < _rank_nodes[r].size(); i++) {
        const CCDGNodeInfo &nd = _rank_nodes[r][i];
        if (nd.type != "SEND" && nd.type != "ISEND" &&
            nd.type != "ALLREDUCE" && nd.type != "BCAST" &&
            nd.type != "GATHER" && nd.type != "ALLGATHER" &&
            nd.type != "SCATTER" && nd.type != "ALLTOALL" &&
            nd.type != "REDUCE") continue;
        uint64_t b = nd.comm_bytes;
        if (b == 0) b = 64; // default minimum message size
        // Collective comm_bytes carries the N-fold total; per packet = /N
        if (nd.type != "SEND" && nd.type != "ISEND") {
          b = b / _num_ranks;
          if (b < 1) b = 1;
        }
        if (b > max_bytes) max_bytes = b;
      }
    }
    _ccdg_inject_queue_depth = _getFlitCount(max_bytes);
    printf("[DBG-INFO] _ccdg_inject_queue_depth : %d", _ccdg_inject_queue_depth);
  }
  cout << "  Injection queue depth: " << _ccdg_inject_queue_depth
       << " flits (PE backpressure)" << endl;

  // Initialize PE states
  _pe_state.resize(_num_ranks);
  for (int r = 0; r < _num_ranks; r++) {
    _pe_state[r].current_node_idx = 0;
    _pe_state[r].remaining_cycles = 0.0;
    _pe_state[r].state = PE_COMPUTE;
    _pe_state[r].pending_edges.clear();
    _pe_state[r].coll_seq = 0;
    _pe_state[r].blocked_cycles = 0.0;
    _pe_state[r].congestion_cycles = 0.0;
    _pe_state[r].compute_cycles_acc = 0.0;
    _pe_state[r].sched_wait_cycles = 0.0;
  }
  _coll_state.resize(_num_ranks);
}

CCDGTrafficManager::~CCDGTrafficManager()
{
}

int CCDGTrafficManager::_getFlitCount(uint64_t bytes) const
{
  int flits = (int)((bytes + _flit_size_bytes - 1) / _flit_size_bytes);
  if (flits < 1) flits = 1;
  return flits;
}

bool CCDGTrafficManager::_parseCCDG(const std::string &filename)
{
  ifstream in(filename.c_str());
  if (!in.is_open()) {
    cerr << "ERROR: Cannot open CCDG file: " << filename << endl;
    return false;
  }

  // Read entire file into string
  string content((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());
  in.close();

  // Simple JSON parser for the CCDG format
  // Expected format:
  // { "num_ranks": N, "nodes": [ {...}, ... ], "cross_rank_edges": [ {...}, ... ] }

  // Find num_ranks
  size_t pos = 0;

  // Parse num_ranks
  pos = content.find("\"num_ranks\"");
  if (pos == string::npos) {
    cerr << "ERROR: num_ranks not found in CCDG file" << endl;
    return false;
  }
  pos = content.find(':', pos);
  _num_ranks = atoi(content.c_str() + pos + 1);
  if (_num_ranks <= 0) {
    cerr << "ERROR: Invalid num_ranks: " << _num_ranks << endl;
    return false;
  }

  // Initialize rank node lists
  _rank_nodes.resize(_num_ranks);

  // Parse nodes array
  pos = content.find("\"nodes\"");
  if (pos == string::npos) {
    cerr << "ERROR: nodes array not found in CCDG file" << endl;
    return false;
  }
  pos = content.find('[', pos);
  if (pos == string::npos) {
    cerr << "ERROR: nodes array opening bracket not found" << endl;
    return false;
  }
  pos++; // skip '['

  // Temporary map from node_id to node index for cross-edge resolution
  map<uint64_t, pair<int, size_t> > node_id_to_idx; // node_id -> (rank, node_index)

  while (pos < content.size()) {
    // Skip whitespace and commas
    while (pos < content.size() && (content[pos] == ' ' || content[pos] == '\n' ||
           content[pos] == '\r' || content[pos] == '\t' || content[pos] == ','))
      pos++;
    if (pos >= content.size() || content[pos] == ']') break;

    if (content[pos] != '{') {
      cerr << "ERROR: Expected '{' at position " << pos << endl;
      return false;
    }
    pos++; // skip '{'

    CCDGNodeInfo node;
    node.id = 0;
    node.rank = 0;
    node.type = "COMPUTE";
    node.compute_cycles = 0.0;
    node.compute_ops = 0.0;
    node.comm_src = -1;
    node.comm_dst = -1;
    node.comm_tag = -1;
    node.comm_bytes = 0;
    node.coll_root = -1;
    node.wse_wavefront_idx = -1;
    node.wse_stage_idx = -1;
    node.wse_kind = "";
    node.wse_mode = "";

    while (pos < content.size()) {
      // Skip whitespace and commas
      while (pos < content.size() && (content[pos] == ' ' || content[pos] == '\n' ||
             content[pos] == '\r' || content[pos] == '\t' || content[pos] == ','))
        pos++;
      if (pos >= content.size() || content[pos] == '}') break;

      // Find key
      if (content[pos] != '"') {
        cerr << "ERROR: Expected key at position " << pos << endl;
        return false;
      }
      pos++; // skip opening quote
      size_t key_end = content.find('"', pos);
      if (key_end == string::npos) {
        cerr << "ERROR: Unterminated key" << endl;
        return false;
      }
      string key = content.substr(pos, key_end - pos);
      pos = key_end + 1;

      // Skip colon
      while (pos < content.size() && content[pos] != ':') pos++;
      pos++; // skip ':'
      while (pos < content.size() && (content[pos] == ' ' || content[pos] == '\n' ||
             content[pos] == '\r' || content[pos] == '\t'))
        pos++;

      // Parse value
      if (key == "id") {
        char *end = NULL;
        node.id = strtoull(content.c_str() + pos, &end, 10);
        pos = end - content.c_str();
      } else if (key == "rank") {
        node.rank = atoi(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else if (key == "type") {
        if (content[pos] == '"') {
          pos++;
          size_t type_end = content.find('"', pos);
          node.type = content.substr(pos, type_end - pos);
          pos = type_end + 1;
        }
      } else if (key == "compute_cycles") {
        node.compute_cycles = atof(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else if (key == "compute_ops") {
        node.compute_ops = atof(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else if (key == "comm_src") {
        node.comm_src = atoi(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else if (key == "comm_dst") {
        node.comm_dst = atoi(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else if (key == "comm_tag") {
        node.comm_tag = atoi(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else if (key == "comm_bytes") {
        char *end = NULL;
        node.comm_bytes = strtoull(content.c_str() + pos, &end, 10);
        pos = end - content.c_str();
      } else if (key == "wse_wavefront_idx") {
        node.wse_wavefront_idx = atoi(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else if (key == "wse_stage_idx") {
        node.wse_stage_idx = atoi(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else if (key == "wse_kind" || key == "wse_mode") {
        if (content[pos] == '"') {
          pos++;
          size_t value_end = content.find('"', pos);
          string value = content.substr(pos, value_end - pos);
          if (key == "wse_kind") node.wse_kind = value;
          else node.wse_mode = value;
          pos = value_end + 1;
        }
      } else if (key == "collective_root") {
        node.coll_root = atoi(content.c_str() + pos);
        while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
               content[pos] != '\n') pos++;
      } else {
        // Skip unknown key's value
        if (content[pos] == '"') {
          pos++;
          size_t val_end = content.find('"', pos);
          if (val_end != string::npos) pos = val_end + 1;
        } else if (content[pos] == '[') {
          int depth = 1;
          pos++;
          while (pos < content.size() && depth > 0) {
            if (content[pos] == '[') depth++;
            else if (content[pos] == ']') depth--;
            pos++;
          }
        } else {
          while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
                 content[pos] != '\n') pos++;
        }
      }
    }
    if (pos < content.size() && content[pos] == '}') pos++; // skip '}'

    // Add node to the appropriate rank
    if (node.rank >= 0 && node.rank < _num_ranks) {
      size_t node_idx = _rank_nodes[node.rank].size();
      _rank_nodes[node.rank].push_back(node);
      node_id_to_idx[node.id] = make_pair(node.rank, node_idx);
    }
  }

  // Parse cross_rank_edges
  pos = content.find("\"cross_rank_edges\"");
  if (pos != string::npos) {
    pos = content.find('[', pos);
    if (pos != string::npos) {
      pos++; // skip '['

      while (pos < content.size()) {
        // Skip whitespace and commas
        while (pos < content.size() && (content[pos] == ' ' || content[pos] == '\n' ||
               content[pos] == '\r' || content[pos] == '\t' || content[pos] == ','))
          pos++;
        if (pos >= content.size() || content[pos] == ']') break;

        if (content[pos] != '{') break;
        pos++; // skip '{'

        CrossRankEdgeInfo edge;
        edge.src_node_id = 0;
        edge.dst_node_id = 0;
        edge.src_rank = -1;
        edge.dst_rank = -1;
        edge.msg_id = -1;
        edge.resolved = false;

        while (pos < content.size()) {
          while (pos < content.size() && (content[pos] == ' ' || content[pos] == '\n' ||
                 content[pos] == '\r' || content[pos] == '\t' || content[pos] == ','))
            pos++;
          if (pos >= content.size() || content[pos] == '}') break;

          if (content[pos] != '"') break;
          pos++;
          size_t key_end = content.find('"', pos);
          string key = content.substr(pos, key_end - pos);
          pos = key_end + 1;

          while (pos < content.size() && content[pos] != ':') pos++;
          pos++;
          while (pos < content.size() && (content[pos] == ' ' || content[pos] == '\n' ||
                 content[pos] == '\r' || content[pos] == '\t'))
            pos++;

          char *end = NULL;
          if (key == "src_node") {
            edge.src_node_id = strtoull(content.c_str() + pos, &end, 10);
            pos = end - content.c_str();
          } else if (key == "dst_node") {
            edge.dst_node_id = strtoull(content.c_str() + pos, &end, 10);
            pos = end - content.c_str();
          } else {
            while (pos < content.size() && content[pos] != ',' && content[pos] != '}' &&
                   content[pos] != '\n') pos++;
          }
        }
        if (pos < content.size() && content[pos] == '}') pos++;

        // Resolve src and dst ranks
        map<uint64_t, pair<int, size_t> >::iterator it;
        it = node_id_to_idx.find(edge.src_node_id);
        if (it != node_id_to_idx.end()) {
          edge.src_rank = it->second.first;
        }
        it = node_id_to_idx.find(edge.dst_node_id);
        if (it != node_id_to_idx.end()) {
          edge.dst_rank = it->second.first;
        }

        size_t edge_idx = _cross_edges.size();
        _cross_edges.push_back(edge);

        // Add edge index to source node's src_edge_indices
        it = node_id_to_idx.find(edge.src_node_id);
        if (it != node_id_to_idx.end()) {
          int r = it->second.first;
          size_t ni = it->second.second;
          if (r >= 0 && r < _num_ranks && ni < _rank_nodes[r].size()) {
            _rank_nodes[r][ni].src_edge_indices.push_back(edge_idx);
          }
        }

        // Add edge index to destination node's dep_edge_indices
        it = node_id_to_idx.find(edge.dst_node_id);
        if (it != node_id_to_idx.end()) {
          int r = it->second.first;
          size_t ni = it->second.second;
          if (r >= 0 && r < _num_ranks && ni < _rank_nodes[r].size()) {
            _rank_nodes[r][ni].dep_edge_indices.push_back(edge_idx);
          }
        }
      }
    }
  }

  cout << "Parsed CCDG: " << _num_ranks << " ranks, "
       << _cross_edges.size() << " cross edges" << endl;
  for (int r = 0; r < _num_ranks; r++) {
    cout << "  Rank " << r << ": " << _rank_nodes[r].size() << " nodes" << endl;
  }

  return true;
}

bool CCDGTrafficManager::_loadScheduleFile()
{
  ifstream in(_ccdg_schedule_file.c_str());
  if (!in.is_open()) {
    cerr << "ERROR: Cannot open schedule file: " << _ccdg_schedule_file << endl;
    return false;
  }

  // node id -> (rank, node index)
  map<uint64_t, pair<int, size_t> > id2idx;
  for (int r = 0; r < (int)_rank_nodes.size(); r++) {
    for (size_t i = 0; i < _rank_nodes[r].size(); i++) {
      id2idx[_rank_nodes[r][i].id] = make_pair(r, i);
    }
  }

  string line;
  int loaded = 0, missing = 0;
  while (getline(in, line)) {
    // strip comments and whitespace
    size_t hash = line.find('#');
    if (hash != string::npos) line = line.substr(0, hash);
    size_t beg = line.find_first_not_of(" \t\r\n");
    if (beg == string::npos) continue;
    line = line.substr(beg);
    istringstream iss(line);
    uint64_t node_id;
    double est;
    if (!(iss >> node_id >> est)) continue;
    map<uint64_t, pair<int, size_t> >::iterator it = id2idx.find(node_id);
    if (it != id2idx.end()) {
      _rank_nodes[it->second.first][it->second.second].sched_est = est;
      loaded++;
    } else {
      missing++;
    }
  }
  in.close();
  cout << "  Schedule: loaded " << loaded << " node release times from "
       << _ccdg_schedule_file;
  if (missing > 0) cout << " (" << missing << " ids not in CCDG)";
  cout << endl;
  return loaded > 0;
}

void CCDGTrafficManager::_injectPacket(int src, int dst, uint64_t bytes, int msg_id)
{
  int flit_count = _getFlitCount(bytes);
  int pid = _cur_pid++;
  int subnetwork = RandomInt(_subnets - 1);

  for (int i = 0; i < flit_count; i++) {
    Flit *f = Flit::New();
    f->id = _cur_id++;
    f->pid = pid;
    f->src = src;
    f->dest = dst;
    f->cl = 0;
    f->head = (i == 0);
    f->tail = (i == flit_count - 1);
    f->type = Flit::ANY_TYPE;
    f->subnetwork = subnetwork;
    f->vc = -1;
    f->pri = 0;
    f->hops = 0;
    f->watch = false;
    f->record = true;
    f->intm = -1;
    f->ph = 0;

    if (i == flit_count - 1) {
      f->data = new uint64_t(msg_id);
    } else {
      f->data = NULL;
    }

    f->ctime = _time;
    f->itime = _time;
    f->atime = _time;

    _total_in_flight_flits[0].insert(make_pair(f->id, f));
    _measured_in_flight_flits[0].insert(make_pair(f->id, f));
    _partial_packets[src][0].push_back(f);
  }
  _total_packets_sent++;
}

void CCDGTrafficManager::_OnPacketIssued(int msg_id, const CCDGNodeInfo &node)
{
  (void)msg_id;
  (void)node;
}

void CCDGTrafficManager::_OnPacketRetired(int msg_id, int dest)
{
  (void)msg_id;
  (void)dest;
}

void CCDGTrafficManager::_collectivePattern(const string &type, int rank, int root,
                                            vector<int> &dests, int &expected_in) const
{
  int N = _num_ranks;
  dests.clear();
  expected_in = 0;

  if (type == "ALLREDUCE" || type == "ALLGATHER") {
    // Recursive doubling; for non power-of-two N use the MPICH-style
    // pre/post steps around the largest power-of-two subset
    int p = 1;
    while (p * 2 <= N) p *= 2;
    if (p == N) {
      for (int mask = 1; mask < N; mask <<= 1) {
        int b = rank ^ mask;
        dests.push_back(b);
        expected_in++;
      }
    } else {
      if (rank >= p) {
        // remainder ranks fold onto rank-p (pre-step send, post-step recv)
        dests.push_back(rank - p);
        expected_in++;
      } else {
        for (int mask = 1; mask < p; mask <<= 1) {
          int b = rank ^ mask;
          dests.push_back(b);
          expected_in++;
        }
        if (rank < N - p) {
          // paired with remainder rank rank+p (pre-step recv, post-step send)
          dests.push_back(rank + p);
          expected_in++;
        }
      }
    }
  } else if (type == "BCAST" || type == "SCATTER") {
    // Binomial tree rooted at root: send to children, recv from parent
    if (root < 0 || root >= N) root = 0;
    int pos = (rank - root + N) % N;
    for (int c = 0; c < 2; c++) {
      int cpos = 2 * pos + 1 + c;
      if (cpos < N) dests.push_back((root + cpos) % N);
    }
    if (pos > 0) expected_in = 1;
  } else if (type == "REDUCE" || type == "GATHER") {
    // Binomial tree reversed: send to parent, recv from children
    if (root < 0 || root >= N) root = 0;
    int pos = (rank - root + N) % N;
    if (pos > 0) dests.push_back((root + (pos - 1) / 2) % N);
    for (int c = 0; c < 2; c++) {
      if (2 * pos + 1 + c < N) expected_in++;
    }
  } else if (type == "ALLTOALL") {
    for (int d = 0; d < N; d++) {
      if (d != rank) {
        dests.push_back(d);
        expected_in++;
      }
    }
  }
  // BARRIER and unknown types: no packets
}

void CCDGTrafficManager::_issueCCDGPacket(int rank, const CCDGNodeInfo &node)
{
  // Only point-to-point sends are handled here; collectives are expanded
  // directly in _advancePE
  if (node.type != "SEND" && node.type != "ISEND") {
    return;
  }

  uint64_t bytes = node.comm_bytes;
  if (bytes == 0 && node.wse_kind != "command") {
    bytes = 64; // default minimum message size
  }

  int dest = node.comm_dst;
  if (dest < 0 || dest >= _nodes) {
    dest = (rank + 1) % _nodes; // fallback
  }

  int msg_id = _msg_id_counter++;

  // Assign msg_id to the cross edges originating from this node
  for (size_t i = 0; i < node.src_edge_indices.size(); i++) {
    size_t eidx = node.src_edge_indices[i];
    if (eidx < _cross_edges.size() && _cross_edges[eidx].msg_id < 0) {
      _cross_edges[eidx].msg_id = msg_id;
      _msg_id_to_edge_idx[msg_id].push_back(eidx);

      // Add dependency to the destination rank
      int dst_rank = _cross_edges[eidx].dst_rank;
      if (dst_rank >= 0 && dst_rank < _num_ranks) {
        _pe_state[dst_rank].pending_edges.insert(eidx);
      }
    }
  }

  _OnPacketIssued(msg_id, node);
  _injectPacket(rank, dest, bytes, msg_id);
}

void CCDGTrafficManager::_advancePE(int rank)
{
  PERankState &pe = _pe_state[rank];
  const vector<CCDGNodeInfo> &nodes = _rank_nodes[rank];

  while (pe.current_node_idx < nodes.size()) {
    const CCDGNodeInfo &node = nodes[pe.current_node_idx];

    // Deterministic orchestration: hold the PE at this node until its
    // compile-time EST (release time). Checked INSIDE the state machine so
    // the dependency wake-up path (_RetireFlit -> _advancePE) cannot bypass
    // release times either - every entry re-checks the current node. A
    // running COMPUTE (remaining cycles being decremented by the main loop)
    // is exempt: it already started inside its window.
    if (_sched_enabled &&
        !(pe.state == PE_COMPUTE && pe.remaining_cycles > 0.0) &&
        node.sched_est > (double)_time) {
      // Park in PE_GATED: the main loop re-enters _advancePE from PE_GATED
      // every cycle, so a release-time wait can never be bypassed. A bare
      // return would leave the PE in PE_COMPUTE with remaining_cycles <= 0,
      // and the main loop would then step current_node_idx once per cycle,
      // skipping every node without injecting a single packet (sent=0).
      pe.sched_wait_cycles += 1.0;
      pe.state = PE_GATED;
      return; // held at this node; caller retries next cycle
    }

    if (node.type == "COMPUTE") {
      // Roofline dwell: compute_ops / compute_rate (ops per NoC cycle).
      // v3 CCDG carries compute_ops (pure-gap CPU time x trace CPU freq);
      // legacy CCDG falls back to compute_cycles / freq_ratio.
      if (node.compute_ops > 0.0) {
        pe.remaining_cycles = node.compute_ops / _ccdg_compute_rate;
      } else {
        // Scale compute cycles by frequency ratio (CPU cycles -> NoC cycles)
        pe.remaining_cycles = node.compute_cycles / _freq_ratio;
      }
      if (pe.remaining_cycles < 1.0) pe.remaining_cycles = 1.0;
      pe.state = PE_COMPUTE;
      return; // We'll process this compute node
    }
    else if (node.type == "SEND" || node.type == "ISEND") {
      // Wavelet gating: hold the PE at this send until the source router
      // is HEAD (and, with VC split, until the stage matches the message
      // direction). The PE retries on every cycle from PE_GATED.
      if (_wse_gating) {
        if (_wse_vc_split && !_wseDirectionOk(rank, node.comm_dst)) {
          pe.state = PE_GATED;
          return;
        }
        if (!_wseCanInject(rank)) {
          pe.state = PE_GATED;
          return;
        }
      }
      // Issue a packet for this send; backpressure first: hold the PE at
      // this send when the injection queue cannot hold the next whole
      // packet (real-hardware credit backpressure stalls the instruction
      // stream; without it congestion is invisible to the CCDG schedule).
      uint64_t queue_bytes =
        (node.comm_bytes || node.wse_kind == "command") ? node.comm_bytes : 64;
      int need = _getFlitCount(queue_bytes);
      if ((int)_partial_packets[rank][0].size() + need > _ccdg_inject_queue_depth) {
        _bp_check_hits++;
        _bp_max_seen = max(_bp_max_seen, (int)_partial_packets[rank][0].size());
        pe.state = PE_BACKPRESSURE;
        return;
      }
      _issueCCDGPacket(rank, node);
      // Move to the next node
      pe.current_node_idx++;
      // Continue the loop to process the next node
    }
    else if (node.type == "RECV" || node.type == "IRECV") {
      // Check if the dependencies for this node are resolved
      bool deps_resolved = true;
      for (size_t i = 0; i < node.dep_edge_indices.size(); i++) {
        size_t eidx = node.dep_edge_indices[i];
        if (eidx < _cross_edges.size() && !_cross_edges[eidx].resolved) {
          deps_resolved = false;
          break;
        }
      }

      if (!deps_resolved) {
        // Block until dependencies are resolved
        pe.state = PE_BLOCKED;
        return;
      }

      // Dependencies resolved, advance to next node
      pe.current_node_idx++;
      // Continue the loop
    }
    else if (node.type == "WAIT" || node.type == "WAITALL" || node.type == "WAITANY") {
      // Check if dependencies (cross-rank edges) are resolved
      bool deps_resolved = true;
      for (size_t i = 0; i < node.dep_edge_indices.size(); i++) {
        size_t eidx = node.dep_edge_indices[i];
        if (eidx < _cross_edges.size() && !_cross_edges[eidx].resolved) {
          deps_resolved = false;
          break;
        }
      }

      if (!deps_resolved) {
        pe.state = PE_BLOCKED;
        return;
      }

      // All dependencies resolved, advance
      pe.current_node_idx++;
      // Continue the loop
    }
    else if (node.type == "ALLREDUCE" || node.type == "BARRIER" ||
             node.type == "BCAST" || node.type == "GATHER" ||
             node.type == "ALLGATHER" || node.type == "SCATTER" ||
             node.type == "ALLTOALL" || node.type == "REDUCE") {
      // Collective operation.
      // BARRIER is a pure global sync point: block until all ranks arrive.
      // Packet-based collectives use a non-blocking traffic model: the
      // expanded packets are injected and the PE advances immediately.
      // (Blocking on collective completion would deadlock whenever the
      // trace-derived wait/send matching is skewed by per-rank wall-clock
      // drift; the CCDG DAG guarantees every SEND eventually fires, and
      // the collective cost still shows up via network contention and the
      // drain of in-flight flits before simulation completion.)
      PECollState &cs = _coll_state[rank];

      if (node.type == "BARRIER") {
        int seq = pe.coll_seq++;
        _barrier_arrived[seq].insert(rank);
        if ((int)_barrier_arrived[seq].size() >= _num_ranks) {
          set<int> arrived = _barrier_arrived[seq];
          _barrier_arrived.erase(seq);
          for (set<int>::iterator it = arrived.begin(); it != arrived.end(); ++it) {
            int r2 = *it;
            _pe_state[r2].current_node_idx++;
            _advancePE(r2);
          }
        } else {
          pe.state = PE_BLOCKED;
        }
        return;
      }

      // Packet-based collective: compute expansion pattern and inject
      // (non-blocking traffic model; per-packet backpressure with partial
      // progress so a stalled collective resumes where it left off).
      vector<int> dests;
      int expected_in = 0;
      _collectivePattern(node.type, rank, node.coll_root, dests, expected_in);
      (void)expected_in;

      // Wavelet gating: hold the PE until the source router is HEAD
      // (no per-destination direction check: collective traffic is
      // omnidirectional by construction)
      if (_wse_gating && !_wseCanInject(rank)) {
        pe.state = PE_GATED;
        return;
      }

      // NOTE: coll_seq is advanced only on first entry (cs.injected == 0);
      // gated/backpressured retries must not consume extra sequence numbers
      // or BARRIER matching desynchronizes across ranks (deadlock).
      if (cs.injected == 0) pe.coll_seq++;

      // Per-packet bytes: comm_bytes carries the N-fold total, divide by N
      uint64_t bytes = node.comm_bytes / _num_ranks;
      if (bytes < 1) bytes = 1;
      int need = _getFlitCount(bytes);

      // Inject expanded packets one at a time; backpressure stalls the PE
      // when the queue cannot hold the next packet. Partial progress is
      // kept in cs.injected and resumed on the next retry.
      for (size_t i = (size_t)cs.injected; i < dests.size(); i++) {
        int dest = dests[i];
        if (dest < 0 || dest >= _nodes || dest == rank) { cs.injected++; continue; }
        if ((int)_partial_packets[rank][0].size() + need > _ccdg_inject_queue_depth) {
          pe.state = PE_BACKPRESSURE;
          return;
        }
        int msg_id = _msg_id_counter++;
        CollMsgInfo ci;
        ci.src = rank;
        ci.dst = dest;
        _coll_msg_info[msg_id] = ci;
        _injectPacket(rank, dest, bytes, msg_id);
        cs.pending_out++;
        cs.injected++;
      }

      // All expanded packets injected: advance to the next node
      cs.injected = 0;
      pe.current_node_idx++;
      // Continue the loop
    }
    else {
      // Unknown node type, skip
      pe.current_node_idx++;
    }
  }

  // All nodes processed
  pe.state = PE_DONE;
  pe.remaining_cycles = 0.0;
}

void CCDGTrafficManager::_RetireFlit(Flit *f, int dest)
{
  // Extract msg_id from tail flit BEFORE calling base class (which frees the flit)
  uint64_t msg_id = 0;
  bool has_msg = (f->tail && f->data != NULL);
  if (has_msg) {
    msg_id = *static_cast<uint64_t *>(f->data);
  }

  // Call base class to handle statistics and free the flit
  TrafficManager::_RetireFlit(f, dest);

  // Resolve dependency if this was a tail flit with a msg_id
  if (has_msg) {
    _OnPacketRetired((int)msg_id, dest);
    // Collective packet: accounting only (non-blocking model)
    map<int, CollMsgInfo>::iterator cit = _coll_msg_info.find((int)msg_id);
    if (cit != _coll_msg_info.end()) {
      CollMsgInfo ci = cit->second;
      _coll_msg_info.erase(cit);
      if (_coll_state[ci.src].pending_out > 0) {
        _coll_state[ci.src].pending_out--;
      }
    } else {
      map<int, vector<size_t> >::iterator it = _msg_id_to_edge_idx.find((int)msg_id);
      if (it != _msg_id_to_edge_idx.end()) {
        // Resolve all edges associated with this msg_id
        for (size_t ei = 0; ei < it->second.size(); ei++) {
          size_t eidx = it->second[ei];
          if (eidx < _cross_edges.size() && !_cross_edges[eidx].resolved) {
            _cross_edges[eidx].resolved = true;

            // Remove from destination rank's pending edges
            int dst_rank = _cross_edges[eidx].dst_rank;
            if (dst_rank >= 0 && dst_rank < _num_ranks) {
              _pe_state[dst_rank].pending_edges.erase(eidx);

              // If the PE was blocked, try to advance it
              if (_pe_state[dst_rank].state == PE_BLOCKED) {
                _advancePE(dst_rank);
              }
            }
          }
        }
      } else {
        static int warn_count = 0;
        if (warn_count++ < 10) {
          cout << "WARNING: msg_id " << msg_id << " not found in map (size="
               << _msg_id_to_edge_idx.size() << ") at time=" << _time
               << " sim_state=" << _sim_state << endl;
        }
      }
    }
  }

  if (has_msg) {
    _total_packets_received++;
  }
}

bool CCDGTrafficManager::_SingleSim()
{
  // Initialize PE states
  for (int r = 0; r < _num_ranks; r++) {
    _pe_state[r].current_node_idx = 0;
    _pe_state[r].pending_edges.clear();
    _pe_state[r].state = PE_COMPUTE;
    _pe_state[r].remaining_cycles = 0.0;
    _pe_state[r].coll_seq = 0;
    _pe_state[r].blocked_cycles = 0.0;
    _pe_state[r].congestion_cycles = 0.0;
    _pe_state[r].compute_cycles_acc = 0.0;
    _coll_state[r] = PECollState();
  }

  // Reset cross edge states
  for (size_t i = 0; i < _cross_edges.size(); i++) {
    _cross_edges[i].msg_id = -1;
    _cross_edges[i].resolved = false;
  }
  _msg_id_to_edge_idx.clear();
  _coll_msg_info.clear();
  _barrier_arrived.clear();
  _msg_id_counter = 0;

  _sim_state = running;
  _time = 0;
  _total_sim_cycles = 0;
  _total_packets_sent = 0;
  _total_packets_received = 0;

  // Initialize all PEs to their first node
  for (int r = 0; r < _num_ranks; r++) {
    _advancePE(r);
  }

  // MAX_CYCLES raised to 2.0e9 (int32-safe) for noc=1GHz runs where
  // compute maps 1:1 to wall time (1 cycle = 1 ns)
  const int MAX_CYCLES = 2000000000;

  while (_time < MAX_CYCLES) {

    // ============ Wavelet (WSE) phase advance ============

    _wseTick();

    // ============ CCDG State Machine ============

    for (int r = 0; r < _num_ranks; r++) {
      PERankState &pe = _pe_state[r];

      if (pe.state == PE_DONE) continue;

      if (pe.state == PE_GATED) {
        // Waiting for an injection window: retry the send/collective
        // (the gating checks inside _advancePE decide whether to proceed)
        _advancePE(r);
        continue;
      }

      if (pe.state == PE_BACKPRESSURE) {
        // Injection queue was full: retry the send/collective. The stall
        // is exposed waiting, tracked separately as congestion_cycles
        // (backpressure) so blocked_cycles stays pure dependency wait.
        pe.congestion_cycles += 1.0;
        _bp_state_cycles++;
        _advancePE(r);
        continue;
      }

      if (pe.state == PE_COMPUTE) {
        // Decrement remaining cycles
        pe.remaining_cycles -= 1.0;
        pe.compute_cycles_acc += 1.0;

        if (pe.remaining_cycles <= 0.0) {
          // Move to the next node
          pe.current_node_idx++;
          _advancePE(r);
        }
      }
      else if (pe.state == PE_BLOCKED) {
        // Waiting on unresolved cross-rank deps (exposed communication).
        // Plan wait outranks dependency wait: while the EST has not been
        // reached the cycle is attributed to the orchestration slot even
        // though the deps are still unresolved (the schedule may be
        // optimistic; the dependency check inside _advancePE still gates
        // correctness).
        const CCDGNodeInfo &n0 = _rank_nodes[r][pe.current_node_idx];
        if (_sched_enabled && n0.sched_est > (double)_time) {
          pe.sched_wait_cycles += 1.0;
        } else {
          pe.blocked_cycles += 1.0;
        }
      }
      // GATED state: do nothing, wait for an injection window
    }

    // ============ Pre-generate packets for CCDG sends ============

    // Note: We handle packet generation in _advancePE, which is called
    // above. The generated packets are in _partial_packets.

    // ============ Network Step ============

    // Read flits and credits from network (same as base _Step)
    bool flits_in_flight = false;
    for (int c = 0; c < _classes; c++) {
      flits_in_flight |= !_total_in_flight_flits[c].empty();
    }
    if (flits_in_flight && (_deadlock_timer++ >= _deadlock_warn_timeout)) {
      _deadlock_timer = 0;
      cout << "WARNING: Possible network deadlock.\n";
    }

    vector<map<int, Flit *> > flits(_subnets);

    for (int subnet = 0; subnet < _subnets; ++subnet) {
      for (int n = 0; n < _nodes; ++n) {
        Flit *const f = _net[subnet]->ReadFlit(n);
        if (f) {
          if (f->watch) {
            *gWatchOut << GetSimTime() << " | "
                       << "node" << n << " | "
                       << "Ejecting flit " << f->id
                       << " (packet " << f->pid << ")"
                       << " from VC " << f->vc
                       << "." << endl;
          }
          flits[subnet].insert(make_pair(n, f));
          if ((_sim_state == warming_up) || (_sim_state == running)) {
            ++_accepted_flits[f->cl][n];
            if (f->tail) {
              ++_accepted_packets[f->cl][n];
            }
          }
        }

        Credit *const c = _net[subnet]->ReadCredit(n);
        if (c) {
#ifdef TRACK_FLOWS
          for (set<int>::const_iterator iter = c->vc.begin(); iter != c->vc.end(); ++iter) {
            int const vc = *iter;
            assert(!_outstanding_classes[n][subnet][vc].empty());
            int cl = _outstanding_classes[n][subnet][vc].front();
            _outstanding_classes[n][subnet][vc].pop();
            assert(_outstanding_credits[cl][subnet][n] > 0);
            --_outstanding_credits[cl][subnet][n];
          }
#endif
          _buf_states[n][subnet]->ProcessCredit(c);
          c->Free();
        }
      }
      _net[subnet]->ReadInputs();
    }

    // Inject packets from _partial_packets into the network
    for (int subnet = 0; subnet < _subnets; ++subnet) {
      for (int n = 0; n < _nodes; ++n) {
        // Wavelet gating: hold flit injection for this node (flit-level
        // window pause; injection resumes on the next HEAD window). Only
        // count cycles where the node actually has queued flits waiting,
        // matching the base _Inject() accounting semantics.
        if (_wse_gating && !_wseCanInject(n)) {
          bool has_pending = false;
          for (int c = 0; c < _classes; ++c) {
            if (!_partial_packets[n][c].empty()) { has_pending = true; break; }
          }
          if (has_pending) ++_wse_inject_blocked_cycles;
          continue;
        }
        Flit *f = NULL;
        BufferState *const dest_buf = _buf_states[n][subnet];
        int const last_class = _last_class[n][subnet];
        int class_limit = _classes;

        if (_hold_switch_for_packet) {
          list<Flit *> const &pp = _partial_packets[n][last_class];
          if (!pp.empty() && !pp.front()->head &&
              !dest_buf->IsFullFor(pp.front()->vc)) {
            f = pp.front();
            assert(f->vc == _last_vc[n][subnet][last_class]);
            --class_limit;
          }
        }

        for (int i = 1; i <= class_limit; ++i) {
          int const c = (last_class + i) % _classes;
          list<Flit *> const &pp = _partial_packets[n][c];

          if (pp.empty()) continue;

          Flit *const cf = pp.front();
          assert(cf);
          assert(cf->cl == c);

          if (cf->subnetwork != subnet) continue;

          if (f && (f->pri >= cf->pri)) continue;

          if (cf->head && cf->vc == -1) {
            OutputSet route_set;
            _rf(NULL, cf, -1, &route_set, true);
            set<OutputSet::sSetElement> const &os = route_set.GetSet();
            assert(os.size() == 1);
            OutputSet::sSetElement const &se = *os.begin();
            assert(se.output_port == -1);
            int vc_start = se.vc_start;
            int vc_end = se.vc_end;
            int vc_count = vc_end - vc_start + 1;
            if (_wse_gating && _wse_vc_split) {
              // Wavelet VC split: restrict to the direction's VC set
              // (4 sets: +x/-x/+y/-y, diagonal round-robins to +/-x)
              int dir = _wseDirClass(n, cf->dest);
              if (dir == 4) dir = (cf->pid % 2) ? 0 : 2;
              int vc_span = vc_end - vc_start + 1;
              int seg = vc_span / 4;
              if (seg >= 1) {
                vc_start = vc_start + dir * seg;
                vc_end = vc_start + seg - 1;
                vc_count = vc_end - vc_start + 1;
              }
            }
            if (_noq) {
              assert(_lookahead_routing);
              const FlitChannel *inject = _net[subnet]->GetInject(n);
              const Router *router = inject->GetSink();
              assert(router);
              int in_channel = inject->GetSinkPort();
              cf->vc = vc_start;
              _rf(router, cf, in_channel, &cf->la_route_set, false);
              cf->vc = -1;
              set<OutputSet::sSetElement> const sl = cf->la_route_set.GetSet();
              assert(sl.size() == 1);
              int next_output = sl.begin()->output_port;
              vc_count /= router->NumOutputs();
              vc_start += next_output * vc_count;
              vc_end = vc_start + vc_count - 1;
            }
            for (int j = 1; j <= vc_count; ++j) {
              int const lvc = _last_vc[n][subnet][c];
              int const vc = (lvc < vc_start || lvc > vc_end) ?
                vc_start : (vc_start + (lvc - vc_start + j) % vc_count);
              assert((vc >= vc_start) && (vc <= vc_end));
              if (!dest_buf->IsAvailableFor(vc)) continue;
              if (dest_buf->IsFullFor(vc)) continue;
              cf->vc = vc;
              break;
            }
          }

          if (cf->vc == -1) continue;
          if (dest_buf->IsFullFor(cf->vc)) continue;
          f = cf;
        }

        if (f) {
          assert(f->subnetwork == subnet);
          int const c = f->cl;

          if (f->head) {
            if (_lookahead_routing) {
              if (!_noq) {
                const FlitChannel *inject = _net[subnet]->GetInject(n);
                const Router *router = inject->GetSink();
                assert(router);
                int in_channel = inject->GetSinkPort();
                _rf(router, f, in_channel, &f->la_route_set, false);
              }
            } else {
              f->la_route_set.Clear();
            }
            dest_buf->TakeBuffer(f->vc);
            _last_vc[n][subnet][c] = f->vc;
          }

          _last_class[n][subnet] = c;
          _partial_packets[n][c].pop_front();

#ifdef TRACK_FLOWS
          ++_outstanding_credits[c][subnet][n];
          _outstanding_classes[n][subnet][f->vc].push(c);
#endif

          dest_buf->SendingFlit(f);

          if (_pri_type == network_age_based) {
            f->pri = numeric_limits<int>::max() - _time;
          }

          if (f->watch) {
            *gWatchOut << GetSimTime() << " | "
                       << "node" << n << " | "
                       << "Injecting flit " << f->id
                       << " into subnet " << subnet
                       << " at time " << _time
                       << " with priority " << f->pri
                       << "." << endl;
          }
          f->itime = _time;

          if (!_partial_packets[n][c].empty() && !f->tail) {
            Flit *const nf = _partial_packets[n][c].front();
            nf->vc = f->vc;
          }

          if ((_sim_state == warming_up) || (_sim_state == running)) {
            ++_sent_flits[c][n];
            if (f->head) {
              ++_sent_packets[c][n];
            }
          }
#ifdef TRACK_FLOWS
          ++_injected_flits[c][n];
#endif
          _net[subnet]->WriteFlit(f, n);
        }
      }
    }

    // Eject flits and write credits
    for (int subnet = 0; subnet < _subnets; ++subnet) {
      for (int n = 0; n < _nodes; ++n) {
        map<int, Flit *>::const_iterator iter = flits[subnet].find(n);
        if (iter != flits[subnet].end()) {
          Flit *const f = iter->second;
          f->atime = _time;
          if (f->watch) {
            *gWatchOut << GetSimTime() << " | "
                       << "node" << n << " | "
                       << "Injecting credit for VC " << f->vc
                       << " into subnet " << subnet
                       << "." << endl;
          }
          Credit *const c = Credit::New();
          c->vc.insert(f->vc);
          _net[subnet]->WriteCredit(c, n);
#ifdef TRACK_FLOWS
          ++_ejected_flits[f->cl][n];
#endif
          _RetireFlit(f, n);
        }
      }
      flits[subnet].clear();
      _net[subnet]->Evaluate();
      _net[subnet]->WriteOutputs();
    }

    ++_time;
    _total_sim_cycles++;

    // ============ Check for completion ============

    bool all_done = true;
    for (int r = 0; r < _num_ranks; r++) {
      if (_pe_state[r].state != PE_DONE) {
        all_done = false;
        break;
      }
    }

    // All PEs done: lift the wavelet gating so that the remaining queued
    // flits can drain (otherwise injection would stall forever)
    if (all_done) {
      _wse_drain_mode = true;
    }

    // Check if all flits have been received
    bool flits_remaining = false;
    for (int c = 0; c < _classes; c++) {
      if (!_total_in_flight_flits[c].empty()) {
        flits_remaining = true;
        break;
      }
    }

    if (all_done && !flits_remaining) {
      break;
    }

    // ============ Dependency deadlock detection ============
    // If no PE is computing, no flit is in flight and no packet is waiting
    // to be injected, blocked PEs can never be woken -> dump diagnostics.
    // NOTE: PE_GATED is NOT a deadlock: the wavelet state machine rotates
    // every wse_phase_width cycles, so a gated send will eventually get an
    // injection window (constructive liveness).
    if (!all_done && !flits_remaining) {
      bool any_compute = false;
      bool any_gated = false;
      bool any_sched = false;
      bool partial_pending = false;
      for (int r = 0; r < _num_ranks; r++) {
        if (_pe_state[r].state == PE_COMPUTE) { any_compute = true; break; }
        if (_pe_state[r].state == PE_GATED) { any_gated = true; break; }
        // A PE parked on its EST release time is alive: release times are
        // outer-clock (cycle counter), not network progress, so the wait
        // always expires when the cycle reaches the EST.
        if (_sched_enabled && _pe_state[r].current_node_idx <
            _rank_nodes[r].size() &&
            _rank_nodes[r][_pe_state[r].current_node_idx].sched_est >
            (double)_time) {
          any_sched = true;
        }
      }
      if (!any_compute && !any_gated && !any_sched) {
        for (int n = 0; n < _nodes && !partial_pending; n++) {
          for (int c = 0; c < _classes; c++) {
            if (!_partial_packets[n][c].empty()) { partial_pending = true; break; }
          }
        }
      }
      if (!any_compute && !any_gated && !any_sched && !partial_pending) {
        cout << "ERROR: CCDG dependency deadlock detected at cycle " << _time << endl;
        cout << "  wse stage=" << _wse_stage << " phase=" << _wse_phase_idx
             << " timer=" << _wse_phase_timer << endl;
        for (int r = 0; r < _num_ranks; r++) {
          PERankState &pe = _pe_state[r];
          cout << "  Rank " << r << ": state=" << pe.state
               << " rstate=" << (int)_router_state[r]
               << " node_idx=" << pe.current_node_idx
               << "/" << _rank_nodes[r].size();
          if (pe.current_node_idx < _rank_nodes[r].size()) {
            const CCDGNodeInfo &nd = _rank_nodes[r][pe.current_node_idx];
            cout << " blocked_on type=" << nd.type << " id=" << nd.id
                 << " dst=" << nd.comm_dst
                 << " dir=" << _wseDirClass(r, nd.comm_dst)
                 << " dep_edges=" << nd.dep_edge_indices.size();
            for (size_t i = 0; i < nd.dep_edge_indices.size(); i++) {
              size_t eidx = nd.dep_edge_indices[i];
              if (eidx < _cross_edges.size() && !_cross_edges[eidx].resolved) {
                cout << " [unresolved e" << eidx
                     << " from rank " << _cross_edges[eidx].src_rank << "]";
              }
            }
          }
          cout << endl;
        }
        break;
      }
    }

    if (gTrace) {
      cout << "TIME " << _time << endl;
    }
  }

  _sim_state = draining;
  _drain_time = _time;

  cout << "CCDG Simulation completed in " << _total_sim_cycles << " cycles." << endl;
  cout << "Backpressure debug: check_hits=" << _bp_check_hits
       << " state_cycles=" << _bp_state_cycles
       << " max_queue_at_trigger=" << _bp_max_seen << endl;
  cout << "Packets sent: " << _total_packets_sent
       << ", received: " << _total_packets_received << endl;

  // Aggregate PE dwell statistics: compute vs blocked (exposed comm wait)
  // vs congestion (send-side backpressure) vs sched_wait (orchestration
  // slots). Conservation: compute+blocked+congestion+sched_wait ==
  // _num_ranks x makespan (gated cycles are the WSE-mode exception).
  double sum_compute = 0.0, sum_blocked = 0.0, sum_congestion = 0.0;
  double sum_sched = 0.0;
  for (int r = 0; r < _num_ranks; r++) {
    sum_compute += _pe_state[r].compute_cycles_acc;
    sum_blocked += _pe_state[r].blocked_cycles;
    sum_congestion += _pe_state[r].congestion_cycles;
    sum_sched += _pe_state[r].sched_wait_cycles;
  }
  double dwell = sum_compute + sum_blocked + sum_congestion + sum_sched;
  double blocked_ratio = (dwell > 0.0) ? (sum_blocked / dwell) : 0.0;
  double congestion_ratio = (dwell > 0.0) ? (sum_congestion / dwell) : 0.0;
  double sched_ratio = (dwell > 0.0) ? (sum_sched / dwell) : 0.0;
  cout << "PE dwell: compute_cycles=" << (long long)sum_compute
       << " blocked_cycles=" << (long long)sum_blocked
       << " congestion_cycles=" << (long long)sum_congestion
       << " sched_wait_cycles=" << (long long)sum_sched
       << " congestion_ratio=" << congestion_ratio
       << " blocked_ratio=" << blocked_ratio
       << " sched_wait_ratio=" << sched_ratio << endl;

  // Check if all cross edges were resolved
  int unresolved = 0;
  for (size_t i = 0; i < _cross_edges.size(); i++) {
    if (!_cross_edges[i].resolved) unresolved++;
  }
  if (unresolved > 0) {
    cout << "WARNING: " << unresolved << " cross-rank edges remain unresolved." << endl;
  }

  return true;
}

void CCDGTrafficManager::_ClearStats()
{
  TrafficManager::_ClearStats();
}

void CCDGTrafficManager::_UpdateOverallStats()
{
  TrafficManager::_UpdateOverallStats();
}

string CCDGTrafficManager::_OverallStatsCSV(int c) const
{
  ostringstream os;
  os << TrafficManager::_OverallStatsCSV(c);
  return os.str();
}

void CCDGTrafficManager::WriteStats(ostream &os) const
{
  TrafficManager::WriteStats(os);
  os << "total_sim_cycles = " << _total_sim_cycles << ";" << endl;
  os << "total_packets_sent = " << _total_packets_sent << ";" << endl;
  os << "total_packets_received = " << _total_packets_received << ";" << endl;
  double sum_compute = 0.0, sum_blocked = 0.0, sum_congestion = 0.0;
  double sum_sched = 0.0;
  for (int r = 0; r < _num_ranks; r++) {
    sum_compute += _pe_state[r].compute_cycles_acc;
    sum_blocked += _pe_state[r].blocked_cycles;
    sum_congestion += _pe_state[r].congestion_cycles;
    sum_sched += _pe_state[r].sched_wait_cycles;
  }
  double dwell = sum_compute + sum_blocked + sum_congestion + sum_sched;
  os << "compute_cycles = " << (long long)sum_compute << ";" << endl;
  os << "blocked_cycles = " << (long long)sum_blocked << ";" << endl;
  os << "congestion_cycles = " << (long long)sum_congestion << ";" << endl;
  os << "sched_wait_cycles = " << (long long)sum_sched << ";" << endl;
  os << "blocked_ratio = " << ((dwell > 0.0) ? sum_blocked / dwell : 0.0) << ";" << endl;
  os << "congestion_ratio = " << ((dwell > 0.0) ? sum_congestion / dwell : 0.0) << ";" << endl;
  os << "sched_wait_ratio = " << ((dwell > 0.0) ? sum_sched / dwell : 0.0) << ";" << endl;
}

void CCDGTrafficManager::DisplayStats(ostream &os) const
{
  TrafficManager::DisplayStats(os);
  os << "Total simulation cycles = " << _total_sim_cycles << endl;
  os << "Total packets sent = " << _total_packets_sent << endl;
  os << "Total packets received = " << _total_packets_received << endl;
  double sum_congestion = 0.0;
  double sum_compute = 0.0, sum_blocked = 0.0, sum_sched = 0.0;
  for (int r = 0; r < _num_ranks; r++) {
    sum_compute += _pe_state[r].compute_cycles_acc;
    sum_blocked += _pe_state[r].blocked_cycles;
    sum_congestion += _pe_state[r].congestion_cycles;
    sum_sched += _pe_state[r].sched_wait_cycles;
  }
  double dwell = sum_compute + sum_blocked + sum_congestion + sum_sched;
  os << "Compute cycles = " << (long long)sum_compute << endl;
  os << "Blocked cycles = " << (long long)sum_blocked << endl;
  os << "congestion_cycles = " << (long long)sum_congestion << endl;
  os << "sched_wait_cycles = " << (long long)sum_sched << endl;
  os << "congestion_ratio = " << ((dwell > 0.0) ? sum_congestion / dwell : 0.0) << endl;
  os << "Blocked ratio = " << ((dwell > 0.0) ? sum_blocked / dwell : 0.0) << endl;
  os << "sched_wait_ratio = " << ((dwell > 0.0) ? sum_sched / dwell : 0.0) << endl;
}

void CCDGTrafficManager::DisplayOverallStats(ostream &os) const
{
  TrafficManager::DisplayOverallStats(os);
  os << "Total simulation cycles = " << _total_sim_cycles << endl;
  os << "Total packets sent = " << _total_packets_sent << endl;
  os << "Total packets received = " << _total_packets_received << endl;
  double sum_compute = 0.0, sum_blocked = 0.0, sum_congestion = 0.0;
  double sum_sched = 0.0;
  for (int r = 0; r < _num_ranks; r++) {
    sum_compute += _pe_state[r].compute_cycles_acc;
    sum_blocked += _pe_state[r].blocked_cycles;
    sum_congestion += _pe_state[r].congestion_cycles;
    sum_sched += _pe_state[r].sched_wait_cycles;
  }
  double dwell = sum_compute + sum_blocked + sum_congestion + sum_sched;
  os << "Compute cycles = " << (long long)sum_compute << endl;
  os << "Blocked cycles = " << (long long)sum_blocked << endl;
  os << "congestion_cycles = " << (long long)sum_congestion << endl;
  os << "sched_wait_cycles = " << (long long)sum_sched << endl;
  os << "Blocked ratio = " << ((dwell > 0.0) ? sum_blocked / dwell : 0.0) << endl;
  os << "congestion_ratio = " << ((dwell > 0.0) ? sum_congestion / dwell : 0.0) << endl;
  os << "sched_wait_ratio = " << ((dwell > 0.0) ? sum_sched / dwell : 0.0) << endl;
}
