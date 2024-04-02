/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: xarray_for_gridbox_data.hpp
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
 * Helpful structs to write data collected from each Gridbox in parallel
 * to an array in a dataset
 */

#ifndef LIBS_OBSERVERS2_XARRAY_FOR_GRIDBOX_DATA_HPP_
#define LIBS_OBSERVERS2_XARRAY_FOR_GRIDBOX_DATA_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <string_view>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/buffer.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"
#include "zarr2/zarr_array.hpp"

/* struct holding an array in a dataset as well a view and its mirror view
which can be useful when collecting data for 1 variable from 'ngbxs' gridboxes
(in parallel) to then write the the array */
template <typename Store, typename T>
struct XarrayForGridboxData {
  XarrayZarrArray<Store, T> xzarr;                         // array in a dataset
  using viewh_data = Buffer<T>::viewh_buffer;              // type of view for h_data
  using mirrorviewd_data = Buffer<T>::mirrorviewd_buffer;  // mirror view type for d_data
  viewh_data h_data;        // view on host for value of 1 variable from every gridbox
  mirrorviewd_data d_data;  // mirror view of h_data on device

  // Constructor to initialize views and pointer to array in dataset
  XarrayForGridboxData(Dataset<Store> &dataset, const std::string_view name,
                       const std::string_view units, const std::string_view dtype,
                       const double scale_factor, const size_t maxchunk, const size_t ngbxs)
      : xzarr(dataset.template create_array<T>(name, units, dtype, scale_factor,
                                               good2Dchunkshape(maxchunk, ngbxs),
                                               {"time", "gbxindex"})),
        h_data("h_data", ngbxs),
        d_data(Kokkos::create_mirror_view(ExecSpace(), h_data)) {}

  // copy data from device view directly to host and then write to array in dataset
  void write_to_array(Dataset<Store> &dataset) {
    Kokkos::deep_copy(h_data, d_data);
    write_to_array(dataset, xzarr, h_data);
  }

  // call function to write shape of array according to dataset
  void write_arrayshape(Dataset<Store> &dataset) { write_arrayshape(dataset, xzarr); }
};

#endif  // LIBS_OBSERVERS2_XARRAY_FOR_GRIDBOX_DATA_HPP_
