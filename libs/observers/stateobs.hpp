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

#include <concepts>
#include <memory>
#include <iostream>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr/statestorage.hpp"

inline Observer auto
StateObserver(const unsigned int interval,
              FSStore &store,
              const int maxchunk,
              const size_t ngbxs);
/* constructs observer of variables in the state 
of each gridbox with a constant timestep 'interval'
using an instance of the DoStateObs class */

struct StateBuffers
{
  std::vector<double> press;
  std::vector<double> temp;
  std::vector<double> qvap;
  std::vector<double> qcond;

  StateBuffers(const unsigned int buffersize)
      : pressbuffer(buffersize, std::numeric_limits<double>::max()),
        tempbuffer(buffersize, std::numeric_limits<double>::max()),
        qvapbuffer(buffersize, std::numeric_limits<double>::max()),
        qcondbuffer(buffersize, std::numeric_limits<double>::max()) {}

  unsigned int copy2buffers(const ThermoState &state, unsigned int j);
  /* copy press, temp, qvap and qcond data in the state to buffers at index j */

  std::pair<unsigned int, unsigned int>
  writechunks(FSStore &store, unsigned int chunkcount);
  /* write buffer vector into attr's store at chunkcount
  and then replace contents of buffer with numeric limit */

  void writejsons(FSStore &store, const std::string &metadata) const;
  /* write same .zarray metadata to a json file for each
  thermostate array in store alongside distinct .zattrs json files */
};

class DoStateObs
/* observe variables in the state of each
gridbox and write them to repspective arrays
in a store as determined by the 
StateStorage instance */
{
private:
  using store_type = StateStorage<double>;
  std::shared_ptr<store_type> zarr;
  
public:
  DoStateObs(FSStore &store,
             const int maxchunk,
             const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                          "<f8", ngbxs))
}

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
