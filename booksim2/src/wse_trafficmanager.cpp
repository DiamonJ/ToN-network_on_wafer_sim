#include "wse_trafficmanager.hpp"

using namespace std;

WSETrafficManager::WSETrafficManager(
    const Configuration &config, const vector<Network *> &net)
    : CCDGTrafficManager(config, net),
      _commands_expected(0), _commands_delivered(0),
      _branches_expected(0), _branches_delivered(0),
      _multicast_branches_delivered(0), _reduction_branches_delivered(0),
      _branch_bytes_expected(0)
{
  _wse_program_file = config.GetStr("wse_program_file");
  for (size_t rank = 0; rank < _rank_nodes.size(); ++rank) {
    for (size_t i = 0; i < _rank_nodes[rank].size(); ++i) {
      const CCDGNodeInfo &node = _rank_nodes[rank][i];
      if (node.wse_wavefront_idx >= 0) {
        _expected_wavefronts.insert(node.wse_wavefront_idx);
      }
      if (node.wse_kind == "command") {
        ++_commands_expected;
      } else if (node.wse_kind == "branch") {
        ++_branches_expected;
        _branch_bytes_expected += node.comm_bytes;
        ++_stage_branches_expected[node.wse_stage_idx];
      }
    }
  }
  cout << "WSETrafficManager: footprint replay wavefronts="
       << _expected_wavefronts.size()
       << " commands=" << _commands_expected
       << " branches=" << _branches_expected
       << " program=" << _wse_program_file << endl;
}

WSETrafficManager::~WSETrafficManager()
{
}

void WSETrafficManager::_OnPacketIssued(int msg_id,
                                        const CCDGNodeInfo &node)
{
  PacketMeta meta;
  meta.wavefront = node.wse_wavefront_idx;
  meta.stage = node.wse_stage_idx;
  meta.kind = node.wse_kind;
  meta.mode = node.wse_mode;
  _packet_meta[msg_id] = meta;
  if (meta.wavefront >= 0) _injected_wavefronts.insert(meta.wavefront);
}

void WSETrafficManager::_OnPacketRetired(int msg_id, int dest)
{
  (void)dest;
  map<int, PacketMeta>::iterator it = _packet_meta.find(msg_id);
  if (it == _packet_meta.end()) return;
  const PacketMeta meta = it->second;
  _packet_meta.erase(it);
  if (meta.kind == "command") {
    ++_commands_delivered;
  } else if (meta.kind == "branch") {
    ++_branches_delivered;
    ++_stage_branches_delivered[meta.stage];
    _stage_completion_cycle[meta.stage] =
      max(_stage_completion_cycle[meta.stage], _time);
    if (meta.mode == "reduction") ++_reduction_branches_delivered;
    else ++_multicast_branches_delivered;
  }
}

double WSETrafficManager::_WSECongestionRatio() const
{
  if (_num_ranks <= 0 || _total_sim_cycles <= 0) return 0.0;
  double congestion = 0.0;
  for (int rank = 0; rank < _num_ranks; ++rank) {
    congestion += _pe_state[rank].congestion_cycles;
  }
  return congestion / ((double)_num_ranks * (double)_total_sim_cycles);
}

bool WSETrafficManager::_SingleSim()
{
  _injected_wavefronts.clear();
  _packet_meta.clear();
  _stage_branches_delivered.clear();
  _stage_completion_cycle.clear();
  _commands_delivered = 0;
  _branches_delivered = 0;
  _multicast_branches_delivered = 0;
  _reduction_branches_delivered = 0;
  bool result = CCDGTrafficManager::_SingleSim();
  cout << "WSE replay: wavefronts=" << _injected_wavefronts.size()
       << "/" << _expected_wavefronts.size()
       << " commands=" << _commands_delivered << "/" << _commands_expected
       << " branches=" << _branches_delivered << "/" << _branches_expected
       << " congestion_ratio=" << _WSECongestionRatio() << endl;
  return result;
}

void WSETrafficManager::WriteStats(ostream &os) const
{
  CCDGTrafficManager::WriteStats(os);
  os << "wse_wavefronts_expected = " << _expected_wavefronts.size() << ";" << endl;
  os << "wse_wavefronts_injected = " << _injected_wavefronts.size() << ";" << endl;
  os << "wse_commands_expected = " << _commands_expected << ";" << endl;
  os << "wse_commands_delivered = " << _commands_delivered << ";" << endl;
  os << "wse_branches_expected = " << _branches_expected << ";" << endl;
  os << "wse_branches_delivered = " << _branches_delivered << ";" << endl;
  os << "wse_multicast_branches_delivered = "
     << _multicast_branches_delivered << ";" << endl;
  os << "wse_reduction_branches_delivered = "
     << _reduction_branches_delivered << ";" << endl;
  os << "wse_branch_bytes_expected = " << _branch_bytes_expected << ";" << endl;
  os << "wse_congestion_ratio = " << _WSECongestionRatio() << ";" << endl;
  for (map<int, uint64_t>::const_iterator it = _stage_branches_expected.begin();
       it != _stage_branches_expected.end(); ++it) {
    int stage = it->first;
    map<int, uint64_t>::const_iterator delivered =
      _stage_branches_delivered.find(stage);
    map<int, int>::const_iterator completion = _stage_completion_cycle.find(stage);
    os << "wse_stage_" << stage << "_branches_expected = " << it->second << ";" << endl;
    os << "wse_stage_" << stage << "_branches_delivered = "
       << (delivered == _stage_branches_delivered.end() ? 0 : delivered->second)
       << ";" << endl;
    os << "wse_stage_" << stage << "_completion_cycle = "
       << (completion == _stage_completion_cycle.end() ? 0 : completion->second)
       << ";" << endl;
  }
}

void WSETrafficManager::DisplayStats(ostream &os) const
{
  CCDGTrafficManager::DisplayStats(os);
  os << "WSE wavefronts injected = " << _injected_wavefronts.size()
     << " / " << _expected_wavefronts.size() << endl;
  os << "WSE branches delivered = " << _branches_delivered
     << " / " << _branches_expected << endl;
  os << "WSE congestion ratio = " << _WSECongestionRatio() << endl;
}
