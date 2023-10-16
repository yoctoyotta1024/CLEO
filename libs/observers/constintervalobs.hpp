/*
 * ----- CLEO -----
 * File: constintervalobs.hpp
 * Project: observers
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 16th October 2023
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

#ifndef CONSTINTERVALOBS_HPP
#define CONSTINTERVALOBS_HPP

#include <iostream>

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "sdmdomain/gridbox.hpp"

struct ConstIntervalObs
{
private:
  unsigned int interval;

  void observe_gbxs(const unsigned int t_mdl,
               const viewh_constgbx h_gbxs) const;

public:
  ConstIntervalObs(const unsigned int obsstep)
      : interval(obsstep) {}

  unsigned int get_obsstep() const
  {
    return interval;
  }

  bool on_step(const unsigned int t_mdl) const
  {
    return (t_mdl % interval == 0);
  }

  void observe_startstep(const unsigned int t_mdl,
                         const viewh_constgbx h_gbxs) const
  /* observe Gridboxes (on host) at start of timestep */
  {
    if (on_step(t_mdl))
    {
      observe_gbxs(t_mdl, h_gbxs);
    }
  }
};

#endif // CONSTINTERVALOBS_HPP