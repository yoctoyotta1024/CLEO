/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: gbxindexobs.hpp
 * Project: observers
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
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

public:
  GbxindexObserver(FSStore &store, const int maxchunk)
      : zarr(std::make_shared<store_type>(store, maxchunk, "gbxindex",
                                          "<u4", " ", 1))
  {
    zarr->is_name("gbxindex");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const
  /* writes gbxindexes to zarr store
  (only if data has not yet been observed) */
  {
    std::cout << "observer includes GbxindexObserver\n";

    if (zarr->get_ndata() == 0)
    {
      const size_t ngbxs(h_gbxs.extent(0));
      for (size_t ii(0); ii < ngbxs; ++ii)
      {
        zarr->value_to_storage(h_gbxs(ii).get_gbxindex());
      }
    }
  }

  void after_timestepping() const {}

  unsigned int next_obs(const unsigned int t_mdl) const
  {
    return LIMITVALUES::uintmax;
  }

  bool on_step(const unsigned int t_mdl) const
  {
    return false;
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {}

  void at_start_step(const unsigned int t_mdl,
                     const Gridbox &gbx) const {}
};

#endif // GBXINDEXOBS_HPP
