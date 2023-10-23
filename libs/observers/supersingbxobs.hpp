/*
 * ----- CLEO -----
 * File: supersingbxobs.hpp
 * Project: observers
 * Created Date: Monday 23rd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 23rd October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 */


#ifndef SUPERSINGBXOBS_HPP
#define SUPERSINGBXOBS_HPP

#include <concepts>
#include <memory>
#include <iostream>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"

inline Observer auto
SupersInGbxObserver(const unsigned int interval,
                    FSStore &store,
                    const int maxchunk);
/* constructs observer of the attributes of 
all superdroplets in each gridbox with a
constant timestep 'interval' using an instance
of the DoStateObs class */

struct DoSupersInGbxObs
{
  
};


inline Observer auto
SupersInGbxObserver(const unsigned int interval,
                    FSStore &store,
                    const int maxchunk)
/* constructs observer of the attributes of 
all superdroplets in each gridbox with a
constant timestep 'interval' using an instance
of the DoStateObs class */
{
  const auto obs = DoSupersInGbxObs(store, maxchunk);
  return ConstTstepObserver(interval, obs);
}
#endif // SUPERSINGBXOBS_HPP