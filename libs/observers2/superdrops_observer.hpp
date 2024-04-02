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
#include <array>
#include <concepts>
#include <iostream>
#include <memory>
#include <string_view>

#include "../cleoconstants.hpp"
#include "./observers.hpp"
#include "./write_gridbox_to_array.hpp"
#include "zarr2/buffer.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"
#include "zarr2/zarr_array.hpp"

// Operator is functor to perform copy of vvel at the centre of each gridbox to d_data
// in parallel. Note conversion of vvel from double (8 bytes) to single precision (4 bytes
// float) in output
struct XiFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto vvel = static_cast<float>(d_gbxs(ii).state.vvelcentre());
    d_data(ii) = vvel;
  }
};

/* constructs observer which writes mass moments of droplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto SuperdropsObserver(const unsigned int interval, Dataset<Store> &dataset,
                                        const int maxchunk, const size_t ngbxs) {
  const WriteGridboxToArray<Store> auto xiwriter = GenericGbxWriter<Store, uint64_t, XiFunc>(
      dataset, "xi", "", "<u8", 1, maxchunk, ngbxs, XiFunc{});

  const auto obsfunc = WriteGridboxes(ParallelGbxsTeamPolicy{}, dataset, xiwriter);

  return ConstTstepObserver(interval, obsfunc);
}

#endif  // LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_
