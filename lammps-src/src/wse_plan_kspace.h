/* -*- c++ -*- ----------------------------------------------------------
   Source-level WSE plan emitter shared by PPPM, Grid3d, and FFT remap.

   Kspace uses a separate per-rank shard from CommBrick.  This avoids
   multiple buffered FILE streams writing the same file while preserving
   the same JSONL record schema and the no-extra-MPI-call guarantee.
------------------------------------------------------------------------- */

#ifndef LMP_WSE_PLAN_KSPACE_H
#define LMP_WSE_PLAN_KSPACE_H

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <string>

namespace LAMMPS_NS {

class WsePlanKspace {
 public:
  static WsePlanKspace &instance()
  {
    static WsePlanKspace writer;
    return writer;
  }

  void begin(int rank, int num_ranks, bool setup, long long timestep)
  {
    rank_ = rank;
    num_ranks_ = num_ranks;
    setup_ = setup;
    timestep_ = timestep;
    active_ = true;
    open();
  }

  void end() { active_ = false; }

  void message(const char *phase, int dst, int value_count, int datatype_bytes)
  {
    if (!active_ || !fp_ || value_count <= 0 || dst == rank_) return;
    const unsigned long long bytes =
      static_cast<unsigned long long>(value_count) *
      static_cast<unsigned long long>(datatype_bytes);
    fprintf(fp_,
            "{\"kind\":\"message\",\"component\":\"kspace\","
            "\"seq\":%lld,\"scope\":\"%s\",\"phase\":\"%s\","
            "\"timestep\":%lld,\"rank\":%d,\"src\":%d,\"dst\":%d,"
            "\"value_count\":%d,\"datatype_bytes\":%d,\"bytes\":%llu}\n",
            seq_++,setup_ ? "setup" : "run",phase,timestep_,rank_,rank_,dst,
            value_count,datatype_bytes,bytes);
  }

  void collective(const char *phase, const char *operation, int value_count,
                  int datatype_bytes)
  {
    if (!active_ || !fp_ || value_count <= 0) return;
    // Match dumpi2ccdg's collective byte convention: one contribution from
    // every rank in the communicator.
    const unsigned long long bytes =
      static_cast<unsigned long long>(num_ranks_) *
      static_cast<unsigned long long>(value_count) *
      static_cast<unsigned long long>(datatype_bytes);
    fprintf(fp_,
            "{\"kind\":\"collective\",\"component\":\"kspace\","
            "\"seq\":%lld,\"scope\":\"%s\",\"phase\":\"%s\","
            "\"operation\":\"%s\",\"timestep\":%lld,\"rank\":%d,"
            "\"value_count\":%d,\"datatype_bytes\":%d,\"bytes\":%llu}\n",
            seq_++,setup_ ? "setup" : "run",phase,operation,timestep_,rank_,
            value_count,datatype_bytes,bytes);
  }

 private:
  WsePlanKspace() :
    fp_(nullptr), seq_(0), rank_(-1), num_ranks_(0), timestep_(0),
    setup_(true), active_(false), checked_(false)
  {}

  ~WsePlanKspace()
  {
    if (fp_) fclose(fp_);
  }

  WsePlanKspace(const WsePlanKspace &) = delete;
  WsePlanKspace &operator=(const WsePlanKspace &) = delete;

  void open()
  {
    if (checked_) return;
    checked_ = true;
    const char *prefix = getenv("LAMMPS_WSE_PLAN");
    if (!prefix || !prefix[0]) return;

    char suffix[48];
    snprintf(suffix,sizeof(suffix),".kspace.rank%04d.jsonl",rank_);
    const std::string filename = std::string(prefix) + suffix;
    fp_ = fopen(filename.c_str(),"w");
    if (!fp_) return;
    setvbuf(fp_,nullptr,_IOFBF,1U << 20);
    fprintf(fp_,
            "{\"kind\":\"metadata\",\"component\":\"kspace\","
            "\"schema_version\":1,\"rank\":%d,\"num_ranks\":%d}\n",
            rank_,num_ranks_);
  }

  FILE *fp_;
  long long seq_;
  int rank_;
  int num_ranks_;
  long long timestep_;
  bool setup_;
  bool active_;
  bool checked_;
};

}    // namespace LAMMPS_NS

#endif
