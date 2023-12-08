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

#include <Kokkos_Core.hpp>

class RunStats
{
private:
  Kokkos::Timer kokkostimer;
  double time1;
  double time2;

public:
  ~RunStats() { summary(); }

  void before_timestepping()
  /* record stats before timestepping
  e.g. current time */
  {
    time1 = kokkostimer.seconds();
  }

  void after_timestepping()
  /* record stats after timestepping
  e.g. current time */
  {
    time2 = kokkostimer.seconds();
  }

  void summary();
  /* print out summary of runtime
  stats to terminal window */
};

#endif // RUNSTATS_HPP