/*
 * ----- CLEO -----
 * File: runstats.hpp
 * Project: observers
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 8th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for making and ouputting
 * statistics related to runtime performance
 * e.g. of timestepping CLEO
 */


#ifndef RUNSTATS_HPP
#define RUNSTATS_HPP 

#include <ios>
#include <iomanip>
#include <iostream>
#include <memory>

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

struct RunStats
{
private:
  Kokkos::Timer kokkostimer;

public:
  double t0; // time of observer creation
  double t_start;
  double t_end;

  RunStats()
      : t0(0.0), t_start(0.0), t_end(0.0),
        kokkostimer(Kokkos::Timer())
  {
    t0 = kokkostimer.seconds(); 
  }

  double time_elapsed() const
  /* returns tiem elapsed since t0 [s] */
  {
    return kokkostimer.seconds() - t0;
  }
};

class RunStatsObserver
{
private:
  unsigned int interval; // timestep between runtime observations
  std::shared_ptr<RunStats> stats;

  void summary() const;
  /* print out summary of runtime
  stats to terminal window */

public:
  RunStatsObserver(const unsigned int obsstep)
      : interval(obsstep), stats(std::make_shared<RunStats>()){}

  void before_timestepping(const viewh_constgbx h_gbxs) const
  /* record stats before timestepping
  e.g. current time */
  {
    stats->t_start = stats->time_elapsed();
  }

  void after_timestepping() const
  /* record stats after timestepping
  e.g. current time */
  {
    stats->t_end = stats->time_elapsed();
    summary();
  }

  unsigned int next_obs(const unsigned int t_mdl) const
  {
    return ((t_mdl / interval) + 1) * interval;
  }

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const
  {
    if (on_step(t_mdl))
    {
      at_start_step();
    }
  }

  void at_start_step(const unsigned int t_mdl,
                     const Gridbox &gbx) const {}
  
  void at_start_step() const {}
};

#endif // RUNSTATS_HPP