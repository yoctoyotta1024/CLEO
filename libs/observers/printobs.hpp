/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: printobs.hpp
 * Project: observers
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 22nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Struct satisfies observer type and has
 * property that observations have a fixed
 * timestep 'interval' between observations
 */

#ifndef LIBS_OBSERVERS_PRINTOBS_HPP_
#define LIBS_OBSERVERS_PRINTOBS_HPP_

#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

namespace dlc = dimless_constants;

struct PrintObserver {
 private:
  unsigned int interval;                              // timestep between print statements
  std::function<double(unsigned int)> step2realtime;  // function to convert timesteps to real time

  void print_statement(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                       const viewd_constsupers totsupers) const;

 public:
  PrintObserver(const unsigned int obsstep, const std::function<double(unsigned int)> step2realtime)
      : interval(obsstep), step2realtime(step2realtime) {}

  void before_timestepping(const viewh_constgbx h_gbxs) const {
    std::cout << "observer includes PrintObserver\n";
  }

  void after_timestepping() const {}

  unsigned int next_obs(const unsigned int t_mdl) const {
    return ((t_mdl / interval) + 1) * interval;
  }

  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  /* observe Gridboxes (on host) at start of timestep */
  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {
    if (on_step(t_mdl)) {
      print_statement(t_mdl, h_gbxs, totsupers);
    }
  }

  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {}
};

#endif  // LIBS_OBSERVERS_PRINTOBS_HPP_
