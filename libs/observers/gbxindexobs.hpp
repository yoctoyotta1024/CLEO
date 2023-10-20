/*
 * ----- CLEO -----
 * File: gbxindexobs.hpp
 * Project: observers
 * Created Date: Friday 20th October 2023
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
 * Observer to output gbxindx to array in a
 * zarr file system storage
 */

#ifndef GBXINDEXOBS_HPP
#define GBXINDEXOBS_HPP

#include <concepts>
#include <memory>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr/coordstorage.hpp"


class GbxindexObserver
/* Observer which makes 1 observation to
observe the gbxindex of each gridbox and
write it to an array 'zarr' store as
determined by the CoordStorage instance */
{
private:
  using store_type = CoordStorage<unsigned int>;
  std::shared_ptr<store_type> zarr;
  bool already_observed;

public:
  GbxindexObserver(FSStore &store, const int maxchunk)
      : zarr(std::make_shared<store_type>(store, maxchunk, "gbxindex",
                                          "<u4", " ", 1)),
        already_observed(false)
  {
    zarr->is_name("gbxindex");
  }

  unsigned int next_obs(const unsigned int t_mdl) const
  {
    return LIMITVALUES::uintmax;
  }

  bool on_step(const unsigned int t_mdl) const
  {
    return false;
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const {}

  void prepare(const viewh_constgbx h_gbxs) const
  /* writes gbxindexes to zarr store
  (only if dat has not yet been observed) */
  {
    if !(already_observed)
    {
      const size_t ngbxs(h_gbxs.extent(0));
      for (size_t ii(0); ii < ngbxs; ++ii)
      {
        zarr.value_to_storage(h_gridboxes(ii).gbxindex);
      }
      already_observed = true;
    }
  }
};

#endif // GBXINDEXOBS_HPP