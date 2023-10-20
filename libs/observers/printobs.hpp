/*
 * ----- CLEO -----
 * File: printobs.hpp
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

#ifndef PRINTOBS_HPP
#define PRINTOBS_HPP

#include <iostream>

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

struct PrintObs
{
private:
  unsigned int interval;

  void observe_start_step(const unsigned int t_mdl,
                          const viewh_constgbx h_gbxs) const;

public:
  PrintObs(const unsigned int obsstep)
      : interval(obsstep) {}

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
      observe_start_step(t_mdl, h_gbxs);
    }
  }
};

#endif // PRINTOBS_HPP