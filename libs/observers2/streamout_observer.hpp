/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: streamout_observer.hpp
 * Project: observers2
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 2nd April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Struct satisfies observer type and streams out live data to an output device
 * (e.g. computer screen) about the state of gridboxes during every observation
 * at fixed 'interval' timesteps.
 */

#ifndef LIBS_OBSERVERS2_STREAMOUT_OBSERVER_HPP_
#define LIBS_OBSERVERS2_STREAMOUT_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

namespace dlc = dimless_constants;

struct StreamOutObserver {
 private:
  unsigned int interval;                              // timestep between print statements
  std::function<double(unsigned int)> step2realtime;  // function to convert timesteps to real time

  void print_statement(const unsigned int t_mdl, const viewd_constgbx d_gbxs) const;

 public:
  StreamOutObserver(const unsigned int obsstep,
                    const std::function<double(unsigned int)> step2realtime)
      : interval(obsstep), step2realtime(step2realtime) {}

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes StreamOutObserver\n";
  }

  void after_timestepping() const {}

  unsigned int next_obs(const unsigned int t_mdl) const {
    return ((t_mdl / interval) + 1) * interval;
  }

  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  /* observe Gridboxes (copy to host) at start of timestep */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs) const {
    if (on_step(t_mdl)) {
      print_statement(t_mdl, d_gbxs);
    }
  }
};

#endif  // LIBS_OBSERVERS2_STREAMOUT_OBSERVER_HPP_
