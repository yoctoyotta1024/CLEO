/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: write_state.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 28th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions to create GridboxDataWriters which write out state variables from
 * each Gridbox, e.g. to use in StateObserver.
 */

#ifndef LIBS_OBSERVERS2_WRITE_STATE_HPP_
#define LIBS_OBSERVERS2_WRITE_STATE_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./write_gridboxes.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"
#include "zarr2/zarr_array.hpp"

// Operator is functor to perform copy of pressure in each gridbox to d_data in parallel
struct PressureFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<double>::mirrorviewd_buffer d_data) const {
    d_data(ii) = d_gbxs(ii).state.press;
  }
};

// returns GridboxDataWriter which writes the pressure in each gridbox to an array in a dataset in a
// store
template <typename Store>
GridboxDataWriter<Store> auto PressureWriter(Dataset<Store> &dataset, const int maxchunk,
                                             const size_t ngbxs) {
  // create shared pointer to 2-D array in dataset for pressure in each gridbox over time
  std::shared_ptr<XarrayZarrArray<Store, double>> xzarr_ptr =
      std::make_shared<XarrayZarrArray<Store, double>>(dataset.template create_array<double>(
          "press", "hPa", "<f8", dlc::P0 / 100, good2Dchunkshape(maxchunk, ngbxs),
          {"time", "gbxindex"}));

  return GenericGbxWriter<Store, double, PressureFunc>(dataset, PressureFunc{}, xzarr_ptr, ngbxs);
}

#endif  // LIBS_OBSERVERS2_WRITE_STATE_HPP_
