/*
 * ----- CLEO -----
 * File: stateobs.hpp
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
 * Observers to output variables from 
 * gridboxes' state to arrays in a zarr
 * file system storage
 */

#ifndef STATEOBS_HPP
#define STATEOBS_HPP

inline Observer auto
StateObserver(const unsigned int interval,
              FSStore &store,
              const int maxchunk,
              const size_t ngbxs);
/* constructs observer of variables in the state 
of each gridbox with a constant timestep 'interval'
using an instance of the DoStateObs class */



inline Observer auto
StateObserver(const unsigned int interval,
              FSStore &store,
              const int maxchunk,
              const size_t ngbxs)
/* constructs observer of variables in the state 
of each gridbox with a constant timestep 'interval'
using an instance of the DoStateObs class */
{
  const auto obs = DoStateObs(store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

#endif // STATEOBS_HPP
