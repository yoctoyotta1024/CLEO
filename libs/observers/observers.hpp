/*
 * ----- CLEO -----
 * File: observers.hpp
 * Project: observers
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 14th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Observer Concept and related structures for various ways
 * of observing (outputing data from) CLEO.
 * An example of an observer is printing some data
 * from a gridbox's thermostate to the terminal
 */


#ifndef OBSERVERS_HPP
#define OBSERVERS_HPP

#include <iostream>

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "sdmdomain/gridbox.hpp"

struct Observer
{
private:
  unsigned int interval;

public:
  Observer(const unsigned int obsstep)
      : interval(obsstep) {}

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void observe(const unsigned int t_mdl,
               const viewh_constgbx h_gbxs) const
  {
    std::cout << "obs @ t = " << t_mdl << "\n";
  }

  void observe_startstep(const unsigned int t_mdl,
                         const viewh_constgbx h_gbxs) const
  /* observe Gridboxes (on host) at start of timestep */
  {
    if (on_step(t_mdl))
    {
      observe(t_mdl, h_gbxs);
    }
  }

  unsigned int get_obsstep() const
  {
    return interval;
  }
};

#endif // OBSERVERS_HPP