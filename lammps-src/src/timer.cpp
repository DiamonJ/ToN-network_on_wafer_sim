/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "timer.h"

#include "comm.h"
#include "error.h"
#include "fmt/chrono.h"
#include "mpi.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   phase trace helpers:
   CLOCK_MONOTONIC nanoseconds, identical clock source to DUMPI wall time
   (see sst-dumpi/dumpi/common/gettime.c), so markers can be correlated
   with per-call timestamps in the DUMPI trace / CCDG nodes.
------------------------------------------------------------------------- */

static const char *phase_names[] = {
    "Total",  "Pair",   "Bond",     "Kspace", "Neigh",   "Comm",   "Modify",
    "Output", "Sync",   "All",      "Dephase", "Dynamics", "Quench", "NEB",
    "RepComm", "RepOut"};

static uint64_t phase_clock_ns()
{
  struct timespec tspec;
  clock_gettime(CLOCK_MONOTONIC, &tspec);
  return (uint64_t) tspec.tv_sec * 1000000000ull + (uint64_t) tspec.tv_nsec;
}

/* simulated clock (MPI_Wtime) for correlating phase markers with
   SMPI/Paje trace events; falls back to real wall time for plain MPI */

static double phase_sim_time()
{
#if defined(MPI_STUBS)
  struct timespec tspec;
  clock_gettime(CLOCK_MONOTONIC, &tspec);
  return (double) tspec.tv_sec + 1e-9 * (double) tspec.tv_nsec;
#else
  return MPI_Wtime();
#endif
}

/* ---------------------------------------------------------------------- */

void (*Timer::phase_marker_hook)(const char *) = Timer::default_phase_marker;
Timer *Timer::_hook_owner = nullptr;

Timer::Timer(LAMMPS *_lmp) : Pointers(_lmp)
{
  _level = NORMAL;
  _sync = OFF;
  _timeout = -1.0;
  _s_timeout = -1.0;
  _checkfreq = 10;
  _nextcheck = -1;
  _phase_fp = nullptr;
  _active_phase = -1;
  this->_stamp(RESET);

  // enable phase trace if requested via environment (rank 0 only)
  const char *trace_path = getenv("LAMMPS_PHASE_TRACE");
  if (trace_path) enable_phase_trace(trace_path);
}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
  disable_phase_trace();
}

/* ----------------------------------------------------------------------
   open the phase trace file (rank 0 only)
------------------------------------------------------------------------- */

void Timer::enable_phase_trace(const char *path)
{
  disable_phase_trace();
  if (comm->me == 0) {
    _phase_fp = fopen(path, "w");
    if (_phase_fp) {
      fprintf(_phase_fp, "# LAMMPS phase trace for DUMPI correlation\n");
      fprintf(_phase_fp, "# clock: CLOCK_MONOTONIC nanoseconds\n");
      fprintf(_phase_fp,
              "# format: walltime_ns,phase_name[,sim_time_s]\n");
      fprintf(_phase_fp,
              "# sim_time_s = MPI_Wtime at the marker (simulated clock "
              "under SMPI)\n");
      fprintf(_phase_fp,
              "# phase rows mark the END of a timer phase: the interval\n");
      fprintf(_phase_fp,
              "# (prev_ts, ts] belongs to the named phase.\n");
      fprintf(_phase_fp,
              "# KSPACE_* rows are START markers of sub-stages inside the\n");
      fprintf(_phase_fp, "# Kspace phase (segment ends at next marker).\n");
      fflush(_phase_fp);
      _hook_owner = this;
    }
  }
  _active_phase = -1;
}

/* ---------------------------------------------------------------------- */

void Timer::disable_phase_trace()
{
  if (_hook_owner == this) _hook_owner = nullptr;
  if (_phase_fp) {
    fclose(_phase_fp);
    _phase_fp = nullptr;
  }
  _active_phase = -1;
}

/* ----------------------------------------------------------------------
   write a named marker at the START of a sub-stage (e.g. inside
   PPPM::compute) so DUMPI calls can be attributed to kspace sub-phases
------------------------------------------------------------------------- */

void Timer::mark_subphase(const char *name)
{
  if (_phase_fp) {
    fprintf(_phase_fp, "%llu,%s,%.9f\n",
            (unsigned long long) phase_clock_ns(), name, phase_sim_time());
    fflush(_phase_fp);
  }
}

/* ----------------------------------------------------------------------
   hook entry point for C-style code (remap.cpp / fft3d.cpp) that cannot
   hold a LAMMPS pointer; forwards to the phase-trace owning Timer
------------------------------------------------------------------------- */

void Timer::default_phase_marker(const char *name)
{
  if (!_hook_owner) return;
  // under SMPI all ranks share this process and call the hook from their
  // own coroutines; only the rank owning the trace file may write
  int me = 0;
  MPI_Comm_rank(_hook_owner->world, &me);
  if (me == 0) _hook_owner->mark_subphase(name);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
  for (int i = 0; i < NUM_TIMER; i++) {
    cpu_array[i] = 0.0;
    wall_array[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void Timer::_stamp(enum ttype which)
{
  double current_cpu = 0.0, current_wall = 0.0;

  if (_level > NORMAL) current_cpu = platform::cputime();
  current_wall = platform::walltime();

  if ((which > TOTAL) && (which < NUM_TIMER)) {
    const double delta_cpu = current_cpu - previous_cpu;
    const double delta_wall = current_wall - previous_wall;

    cpu_array[which] += delta_cpu;
    wall_array[which] += delta_wall;
    cpu_array[ALL] += delta_cpu;
    wall_array[ALL] += delta_wall;

    // phase trace: record end of this timer phase (rank 0 file only)
    if (_phase_fp && (which != _active_phase)) {
      _active_phase = which;
      fprintf(_phase_fp, "%llu,%s,%.9f\n",
              (unsigned long long) phase_clock_ns(), phase_names[which],
              phase_sim_time());
      fflush(_phase_fp);
    }
  }

  previous_cpu = current_cpu;
  previous_wall = current_wall;

  if (which == RESET) {
    this->init();
    cpu_array[TOTAL] = current_cpu;
    wall_array[TOTAL] = current_wall;
  }

  if (_sync) {
    MPI_Barrier(world);
    if (_level > NORMAL) current_cpu = platform::cputime();
    current_wall = platform::walltime();

    cpu_array[SYNC] += current_cpu - previous_cpu;
    wall_array[SYNC] += current_wall - previous_wall;
    previous_cpu = current_cpu;
    previous_wall = current_wall;
  }
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start()
{
  double current_cpu = 0.0, current_wall = 0.0;

  MPI_Barrier(world);

  if (_level < LOOP) return;

  current_cpu = platform::cputime();
  current_wall = platform::walltime();

  cpu_array[TOTAL] = current_cpu;
  wall_array[TOTAL] = current_wall;
  previous_cpu = current_cpu;
  previous_wall = current_wall;
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop()
{
  double current_cpu = 0.0, current_wall = 0.0;

  MPI_Barrier(world);

  if (_level < LOOP) return;

  current_cpu = platform::cputime();
  current_wall = platform::walltime();

  cpu_array[TOTAL] = current_cpu - cpu_array[TOTAL];
  wall_array[TOTAL] = current_wall - wall_array[TOTAL];
}

/* ---------------------------------------------------------------------- */

double Timer::cpu(enum ttype which)
{
  double current_cpu = platform::cputime();
  return (current_cpu - cpu_array[which]);
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(enum ttype which)
{
  if (_level == OFF) return 0.0;
  double current_wall = platform::walltime();
  return (current_wall - wall_array[which]);
}

/* ---------------------------------------------------------------------- */

void Timer::set_wall(enum ttype which, double newtime)
{
  wall_array[which] = newtime;
}

/* ---------------------------------------------------------------------- */

void Timer::init_timeout()
{
  _s_timeout = _timeout;
  if (_timeout < 0)
    _nextcheck = -1;
  else
    _nextcheck = _checkfreq;
}

/* ---------------------------------------------------------------------- */

void Timer::print_timeout(FILE *fp)
{
  if (!fp) return;

  // format timeout setting
  if (_timeout > 0) {
    // time since init_timeout()
    const double d = platform::walltime() - timeout_start;
    // remaining timeout in seconds
    int s = _timeout - d;
    // remaining 1/100ths of seconds
    const int hs = 100 * ((_timeout - d) - s);
    // breaking s down into second/minutes/hours
    const int seconds = s % 60;
    s = (s - seconds) / 60;
    const int minutes = s % 60;
    const int hours = (s - minutes) / 60;
    fprintf(fp, "  Walltime left : %d:%02d:%02d.%02d\n", hours, minutes, seconds, hs);
  }
}

/* ---------------------------------------------------------------------- */

bool Timer::_check_timeout()
{
  double walltime = platform::walltime() - timeout_start;
  // broadcast time to ensure all ranks act the same.
  MPI_Bcast(&walltime, 1, MPI_DOUBLE, 0, world);

  if (walltime < _timeout) {
    _nextcheck += _checkfreq;
    return false;
  } else {
    if (comm->me == 0) error->warning(FLERR, "Wall time limit reached");
    _timeout = 0.0;
    return true;
  }
}

/* ---------------------------------------------------------------------- */
double Timer::get_timeout_remain()
{
  double remain = _timeout + timeout_start - platform::walltime();
  // never report a negative remaining time.
  if (remain < 0.0) remain = 0.0;
  return (_timeout < 0.0) ? 0.0 : remain;
}

/* ----------------------------------------------------------------------
   modify parameters of the Timer class
------------------------------------------------------------------------- */
static const char *timer_style[] = {"off", "loop", "normal", "full"};
static const char *timer_mode[] = {"nosync", "(dummy)", "sync"};

void Timer::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], timer_style[OFF]) == 0) {
      _level = OFF;
    } else if (strcmp(arg[iarg], timer_style[LOOP]) == 0) {
      _level = LOOP;
    } else if (strcmp(arg[iarg], timer_style[NORMAL]) == 0) {
      _level = NORMAL;
    } else if (strcmp(arg[iarg], timer_style[FULL]) == 0) {
      _level = FULL;
    } else if (strcmp(arg[iarg], timer_mode[OFF]) == 0) {
      _sync = OFF;
    } else if (strcmp(arg[iarg], timer_mode[NORMAL]) == 0) {
      _sync = NORMAL;
    } else if (strcmp(arg[iarg], "timeout") == 0) {
      ++iarg;
      if (iarg < narg) {
        _timeout = utils::timespec2seconds(arg[iarg]);
      } else
        error->all(FLERR, "Illegal timer command");
    } else if (strcmp(arg[iarg], "every") == 0) {
      ++iarg;
      if (iarg < narg) {
        _checkfreq = utils::inumeric(FLERR, arg[iarg], false, lmp);
        if (_checkfreq <= 0) error->all(FLERR, "Illegal timer command");
      } else
        error->all(FLERR, "Illegal timer command");
    } else
      error->all(FLERR, "Illegal timer command");
    ++iarg;
  }

  timeout_start = platform::walltime();
  if (comm->me == 0) {

    // format timeout setting
    std::string timeout = "off";
    if (_timeout >= 0.0) {
      std::tm tv = fmt::gmtime((std::time_t) _timeout);
      timeout = fmt::format("{:02d}:{:%M:%S}", tv.tm_yday * 24 + tv.tm_hour, tv);
    }

    utils::logmesg(lmp, "New timer settings: style={}  mode={}  timeout={}\n", timer_style[_level],
                   timer_mode[_sync], timeout);
  }
}
