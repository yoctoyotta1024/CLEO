/*
 * ----- CLEO -----
 * File: runtimestats.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * class for making and ouputting statistics related to runtime
 * whilst timestepping CLEO
 */


#ifndef RUNSTATS_HPP 
#define RUNSTATS_HPP 

#include <iostream>

#include <Kokkos_Core.hpp>

class RunStats
{
private:
  Kokkos::Timer kokkostimer;
  double time1;
  double time2;

public:
  void pre_timestepping()
  {
    time1 = kokkostimer.seconds();
  }

  void post_timestepping()
  {
    time2 = kokkostimer.seconds();
  }

  void summary();
};

#endif // RUNSTATS_HPP