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
struct PressFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto press = static_cast<float>(d_gbxs(ii).state.press);
    d_data(ii) = press;
  }
};

// returns GridboxDataWriter which writes the pressure in each gridbox to an array in a dataset in a
// store
template <typename Store>
GridboxDataWriter<Store> auto PressWriter(Dataset<Store> &dataset, const int maxchunk,
                                          const size_t ngbxs) {
  // create shared pointer to 2-D array in dataset for pressure in each gridbox over time
  std::shared_ptr<XarrayZarrArray<Store, float>> xzarr_ptr =
      std::make_shared<XarrayZarrArray<Store, float>>(dataset.template create_array<float>(
          "press", "hPa", "<f4", dlc::P0 / 100, good2Dchunkshape(maxchunk, ngbxs),
          {"time", "gbxindex"}));

  return GenericGbxWriter<Store, float, PressFunc>(dataset, PressFunc{}, xzarr_ptr, ngbxs);
}

// Operator is functor to perform copy of temperature in each gridbox to d_data in parallel
struct TempFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto temp = static_cast<float>(d_gbxs(ii).state.temp);
    d_data(ii) = temp
  }
};

// returns GridboxDataWriter which writes the temperature in each gridbox to an array in a
// dataset in a store
template <typename Store>
GridboxDataWriter<Store> auto TempWriter(Dataset<Store> &dataset, const int maxchunk,
                                         const size_t ngbxs) {
  // create shared pointer to 2-D array in dataset for temperature in each gridbox over time
  std::shared_ptr<XarrayZarrArray<Store, float>> xzarr_ptr =
      std::make_shared<XarrayZarrArray<Store, float>>(dataset.template create_array<float>(
          "temp", "K", "<f4", dlc::TEMP0, good2Dchunkshape(maxchunk, ngbxs), {"time", "gbxindex"}));

  return GenericGbxWriter<Store, float, TempFunc>(dataset, TempFunc{}, xzarr_ptr, ngbxs);
}

#endif  // LIBS_OBSERVERS2_WRITE_STATE_HPP_
