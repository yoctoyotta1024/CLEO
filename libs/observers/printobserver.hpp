/*
 * ----- CLEO -----
 * File: printobserver.hpp
 * Project: observers
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Struct satisfies observer type and has
 * property that observations have a fixed
 * timestep 'interval' between observations
 */

#ifndef PRINTOBSERVER_HPP
#define PRINTOBSERVER_HPP

#include <ios>
#include <iostream>
#include <iomanip>
#include <functional>

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "../cleoconstants.hpp"
#include "gridboxes/gridbox.hpp"

namespace dlc = dimless_constants;

struct PrintObserver
{
private:
  unsigned int interval;                          // timestep between print statements
  std::function<double(int)> step2realtime; // function to convert timesteps to real time

  void print_statement(const unsigned int t_mdl,
                          const viewh_constgbx h_gbxs) const;

public:
  PrintObserver(const unsigned int obsstep,
                const std::function<double(int)> step2realtime)
      : interval(obsstep), step2realtime(step2realtime) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int next_obs(const unsigned int t_mdl) const
  {
    return ((t_mdl / interval) + 1) * interval;
  }

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const
  /* observe Gridboxes (on host) at start of timestep */
  {
    if (on_step(t_mdl))
    {
      print_statement(t_mdl, h_gbxs);
    }
  }
};

#endif // PRINTOBSERVER_HPP