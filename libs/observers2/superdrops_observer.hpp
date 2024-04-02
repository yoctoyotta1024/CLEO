/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: superdrops_observer.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 2nd April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output variables related to Gridboxes' state at the start of
 * each timestep to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "./do_write_supers.hpp"
#include "./observers.hpp"
#include "./write_gridbox_to_array.hpp"
#include "zarr2/dataset.hpp"

/* Operator is functor to perform copy of xi for each superdroplet in each gridbox to d_data
in parallel. Note conversion of xi from size_t (arch dependent usuualyl 8 bytes) to long
precision unsigned int (unit64_t) */
struct XiFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamMember &team_member, viewd_constgbx d_gbxs,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    const int ii = team_member.league_rank();
    auto supers(d_gbxs(ii).supersingbx.readonly());
    const size_t kk = 0;
    d_data(ii) = static_cast<uint32_t>(supers(kk).get_xi());  // TODO(CB) WIP
  }
};

/* returns WriteGridboxToArray which writes the number of superdrops in each
gridbox to an array in a dataset in a store called "raggedcount_nsupers" */
template <typename Store>
WriteGridboxToArray<Store> auto RaggedCountNsupersWriter(const Dataset<Store> &dataset,
                                                         const size_t maxchunk,
                                                         const size_t ngbxs) {
  return GenericWriteGridboxToXarray<Store, uint32_t, NsupersFunc>(
      dataset, "raggedcount_nsupers", "", "<u4", 1, maxchunk, ngbxs, NsupersFunc{});
}

/* constructs observer which writes mass moments of droplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto SuperdropsObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                        const int maxchunk, const size_t ngbxs) {
  const WriteGridboxToArray<Store> auto xi = GenericWriteGridboxToXarray<Store, uint32_t, XiFunc>(
      dataset, "xi", "", "<u8", 1, maxchunk, ngbxs, XiFunc{});

  const WriteGridboxToArray<Store> auto raggedcount =
      RaggedCountNsupersWriter(dataset, maxchunk, ngbxs);

  const auto obsfunc =
      DoWriteSupers(ParallelGbxsTeamPolicy{}, dataset, CombineWG2A<Store>{}(raggedcount, xi));
  return ConstTstepObserver(interval, obsfunc);
}

#endif  // LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_
