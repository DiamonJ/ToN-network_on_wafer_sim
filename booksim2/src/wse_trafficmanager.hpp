#ifndef _WSE_TRAFFICMANAGER_HPP_
#define _WSE_TRAFFICMANAGER_HPP_

#include <map>
#include <set>
#include <string>

#include "ccdg_trafficmanager.hpp"

// Phase-1 manager-level multicast replay. The compiler lowers each WSE root
// injection to a synchronization-free scheduled carrier; the proven tree
// footprint remains in the sibling .program.json artifact.
class WSETrafficManager : public CCDGTrafficManager {
protected:
  struct PacketMeta {
    int wavefront;
    int stage;
    std::string kind;
    std::string mode;
  };

  std::set<int> _expected_wavefronts;
  std::set<int> _injected_wavefronts;
  std::map<int, PacketMeta> _packet_meta;
  std::map<int, uint64_t> _stage_branches_expected;
  std::map<int, uint64_t> _stage_branches_delivered;
  std::map<int, int> _stage_completion_cycle;
  std::string _wse_program_file;
  uint64_t _commands_expected;
  uint64_t _commands_delivered;
  uint64_t _branches_expected;
  uint64_t _branches_delivered;
  uint64_t _multicast_branches_delivered;
  uint64_t _reduction_branches_delivered;
  uint64_t _branch_bytes_expected;

  virtual void _OnPacketIssued(int msg_id, const CCDGNodeInfo &node);
  virtual void _OnPacketRetired(int msg_id, int dest);
  virtual bool _SingleSim();
  double _WSECongestionRatio() const;

public:
  WSETrafficManager(const Configuration &config,
                    const std::vector<Network *> &net);
  virtual ~WSETrafficManager();
  virtual void WriteStats(std::ostream &os = std::cout) const;
  virtual void DisplayStats(std::ostream &os = std::cout) const;
};

#endif
