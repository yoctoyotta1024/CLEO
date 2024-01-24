/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: runstats.hpp
 * Project: observers
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 9th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for making and ouputting
 * statistics related to runtime performance
 * e.g. of timestepping CLEO
 */

#ifndef LIBS_OBSERVERS_RUNSTATS_HPP_
#define LIBS_OBSERVERS_RUNSTATS_HPP_

#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

struct RunStats {
 private:
  Kokkos::Timer kokkostimer;

 public:
  double t0;  // time of observer creation
  double t_start;
  double t_end;

  RunStats() : kokkostimer(Kokkos::Timer()), t0(0.0), t_start(0.0), t_end(0.0) {
    t0 = kokkostimer.seconds();
  }

  /* returns tiem elapsed since t0 [s] */
  double time_elapsed() const { return kokkostimer.seconds() - t0; }
};

class RunStatsObserver {
 private:
  unsigned int interval;  // timestep between runtime observations
  std::shared_ptr<RunStats> stats;
  std::string stats_filename;

  /* print out summary of runtime
  stats to terminal window */
  void print_summary() const;

  /* writes out some of the runtime
  statistics to the statsfile */
  void write_to_file() const;

 public:
  RunStatsObserver(const unsigned int obsstep, const std::string stats_filename)
      : interval(obsstep), stats(std::make_shared<RunStats>()), stats_filename(stats_filename) {}

  /* record stats before timestepping
  e.g. current time */
  void before_timestepping(const viewh_constgbx h_gbxs) const {
    stats->t_start = stats->time_elapsed();
  }

  /* record stats after timestepping
  e.g. current time */
  void after_timestepping() const {
    stats->t_end = stats->time_elapsed();
    print_summary();
    write_to_file();
  }

  unsigned int next_obs(const unsigned int t_mdl) const {
    return ((t_mdl / interval) + 1) * interval;
  }

  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {
    if (on_step(t_mdl)) {
      at_start_step();
    }
  }

  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {}

  void at_start_step() const {}
};

#endif  // LIBS_OBSERVERS_RUNSTATS_HPP_
