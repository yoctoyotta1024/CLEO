/*
 * ----- CLEO -----
 * File: supersattrsobs.hpp
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

#ifndef SUPERSATTRSOBS_HPP
#define SUPERSATTRSOBS_HPP

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
SupersAttrsObserver(const unsigned int interval,
                    FSStore &store,
                    const int maxchunk);
/* constructs observer of the attributes of
all superdroplets in each gridbox with a
constant timestep 'interval' using an instance
of the DoStateObs class */

template <typename ContiguousRaggedStorage>
struct DoSupersAttrsObs
/* observe superdroplets by writing their (attributes')
data to contigious ragged represented arrays as
determined by the ContiguousRaggedSDStorage instance */
{
private:
  std::shared_ptr<ContiguousRaggedStorage> zarr;

public:
  DoSupersAttrsObs(FSStore &store,
                   const int maxchunk,
                   const size_t ngbxs)
      : zarr(std::make_shared<ContiguousRaggedStorage>()) {}
  
  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    std::cout << "observer includes SupersAttrsObserver\n";
  } 

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const
  /* writes some variables from gridbox state
  to 2-D zarr storages as determined by the
  StateBuffers struct */
  {
    const size_t ngbxs(h_gbxs.extent(0));
    size_t totnsupers(0);
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      auto h_supers = h_gbxs(ii).supersingbx.hostcopy();

      for (size_t kk(0); kk < h_supers.extent(0); ++kk)
      {
        const auto superdrop = h_supers(kk);
        zarr->data_to_raggedstorage(superdrop);
        ++totnsupers;
      }
    }
    zarr -> raggedarray_count(totnsupers);
  }
};

inline Observer auto
SupersAttrsObserver(const unsigned int interval,
                    FSStore &store,
                    const int maxchunk)
/* constructs observer of the attributes of
all superdroplets in each gridbox with a
constant timestep 'interval' using an instance
of the DoStateObs class */
{
  const auto obs = DoSupersAttrsObs(store, maxchunk);
  return ConstTstepObserver(interval, obs);
}
#endif // SUPERSATTRSOBS_HPP