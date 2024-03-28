/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: state_writers.hpp
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

#ifndef LIBS_OBSERVERS2_STATE_WRITERS_HPP_
#define LIBS_OBSERVERS2_STATE_WRITERS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>
#include <string_view>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./write_gridboxes.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"
#include "zarr2/zarr_array.hpp"

// Operator is functor to perform copy of pressure in each gridbox to d_data in parallel.
// Note conversion of pressure from double (8 bytes) to single precision (4 bytes float) in output
struct PressFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto press = static_cast<float>(d_gbxs(ii).state.press);
    d_data(ii) = press;
  }
};

// Operator is functor to perform copy of temperature in each gridbox to d_data in parallel.
// Note conversion of temperature from double (8 bytes) to single precision (4 bytes float) in
// output
struct TempFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto temp = static_cast<float>(d_gbxs(ii).state.temp);
    d_data(ii) = temp;
  }
};

// Operator is functor to perform copy of vapour mass mixing ratio (qvap) in each gridbox to d_data
// in parallel. Note conversion of qvap from double (8 bytes) to single precision (4 bytes float) in
// output
struct QvapFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto qvap = static_cast<float>(d_gbxs(ii).state.qvap);
    d_data(ii) = qvap;
  }
};

// Operator is functor to perform copy of liquid mass mixing ratio (qcond) in each gridbox to d_data
// in parallel. Note conversion of qcond from double (8 bytes) to single precision (4 bytes
// float) in output
struct QcondFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto qcond = static_cast<float>(d_gbxs(ii).state.qcond);
    d_data(ii) = qcond;
  }
};

// Operator is functor to perform copy of wvel at the centre of each gridbox to d_data
// in parallel. Note conversion of wvel from double (8 bytes) to single precision (4 bytes
// float) in output
struct WvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto wvel = static_cast<float>(d_gbxs(ii).state.wvelcentre());
    d_data(ii) = wvel;
  }
};

// Operator is functor to perform copy of uvel at the centre of each gridbox to d_data
// in parallel. Note conversion of uvel from double (8 bytes) to single precision (4 bytes
// float) in output
struct UvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto uvel = static_cast<float>(d_gbxs(ii).state.uvelcentre());
    d_data(ii) = uvel;
  }
};

// Operator is functor to perform copy of vvel at the centre of each gridbox to d_data
// in parallel. Note conversion of vvel from double (8 bytes) to single precision (4 bytes
// float) in output
struct VvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto vvel = static_cast<float>(d_gbxs(ii).state.vvelcentre());
    d_data(ii) = vvel;
  }
};

// returns GridboxDataWriter which writes the pressure, temperature, qvap, and qcond from
// each gridbox to an array in a dataset in a store
template <typename Store>
GridboxDataWriter<Store> auto ThermoWriter(Dataset<Store> &dataset, const int maxchunk,
                                           const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);

  auto make_array_ptr =
      [&dataset, &chunkshape](
          const std::string_view name, const std::string_view units,
          const double scale_factor) -> std::shared_ptr<XarrayZarrArray<Store, float>> {
    return std::make_shared<XarrayZarrArray<Store, float>>(dataset.template create_array<float>(
        name, units, "<f4", scale_factor, chunkshape, {"time", "gbxindex"}));
  };

  // create shared pointer to 2-D array in dataset for pressure in each gridbox over time
  auto press_ptr = make_array_ptr("press", "hPa", dlc::P0 / 100);

  // create shared pointer to 2-D array in dataset for temperature in each gridbox over time
  auto temp_ptr = make_array_ptr("temp", "K", dlc::TEMP0);

  // create shared pointer to 2-D array in dataset for qvap in each gridbox over time
  auto qvap_ptr = make_array_ptr("qvap", "g/Kg", 1000.0);

  // create shared pointer to 2-D array in dataset for qcond in each gridbox over time
  auto qcond_ptr = make_array_ptr("qcond", "g/Kg", 1000.0);

  const auto c = CombineGDW<Store>{};
  auto press = GenericGbxWriter<Store, float, PressFunc>(dataset, PressFunc{}, press_ptr, ngbxs);
  auto temp = GenericGbxWriter<Store, float, TempFunc>(dataset, TempFunc{}, temp_ptr, ngbxs);
  auto qvap = GenericGbxWriter<Store, float, QvapFunc>(dataset, QvapFunc{}, qvap_ptr, ngbxs);
  auto qcond = GenericGbxWriter<Store, float, QcondFunc>(dataset, QcondFunc{}, qcond_ptr, ngbxs);

  return c(c(qvap, c(press, temp)), qcond);
}

// returns GridboxDataWriter which writes the wind velocity components from the centre of each
// gridbox to an array in a dataset in a store
template <typename Store>
GridboxDataWriter<Store> auto WindVelocityWriter(Dataset<Store> &dataset, const int maxchunk,
                                                 const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);

  auto vel_ptr =
      [&dataset,
       &chunkshape](const std::string_view name) -> std::shared_ptr<XarrayZarrArray<Store, float>> {
    return std::make_shared<XarrayZarrArray<Store, float>>(dataset.template create_array<float>(
        name, "m/s", "<f4", dlc::W0, chunkshape, {"time", "gbxindex"}));
  };

  // create shared pointer to 2-D arrays in a datasetfor the velocity at the centre of each gridbox
  // over time and use to make GbxWriters for each velocity component
  auto wvel = GenericGbxWriter<Store, float, WvelFunc>(dataset, WvelFunc{}, vel_ptr("wvel"), ngbxs);
  auto uvel = GenericGbxWriter<Store, float, UvelFunc>(dataset, UvelFunc{}, vel_ptr("uvel"), ngbxs);
  auto vvel = GenericGbxWriter<Store, float, VvelFunc>(dataset, VvelFunc{}, vel_ptr("vvel"), ngbxs);

  const auto c = CombineGDW<Store>{};
  return c(wvel, c(vvel, uvel));
}

#endif  // LIBS_OBSERVERS2_STATE_WRITERS_HPP_
