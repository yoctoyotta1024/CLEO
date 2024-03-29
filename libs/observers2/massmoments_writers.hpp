/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: massmoments_writers.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 29th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions to create GridboxDataWriters which write out state variables from
 * each Gridbox, e.g. to use in StateObserver.
 */

#ifndef LIBS_OBSERVERS2_MASSMOMENTS_WRITERS_HPP_
#define LIBS_OBSERVERS2_MASSMOMENTS_WRITERS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>
#include <string_view>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./write_gridboxes.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"
#include "zarr2/zarr_array.hpp"

/* Operator is functor to calculate 0th mass moment (i.e. 0th radius moment, i.e. number of
droplets) in each gridbox to d_data in parallel. Note conversion of moment from double (8 bytes) to
single precision (4 bytes float) in output.
*/
struct MassMom0Func {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto m0 = (d_gbxs(ii).state.qcond);
    d_data(ii) = static_cast<float>(m0);  // TODO(CB)
  }
};

/* returns GridboxDataWriter which writes the 0th, 1st and 2nd mass moment
(i.e. 0th, 3rd and 6th radius moment) from each gridbox to arrays in a dataset in a store
*/
template <typename Store>
GridboxDataWriter<Store> auto MassMomentsWriter(Dataset<Store> &dataset, const int maxchunk,
                                                const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);

  auto make_array_ptr =
      [&dataset, &chunkshape](
          const std::string_view name, const std::string_view units,
          const double scale_factor) -> std::shared_ptr<XarrayZarrArray<Store, float>> {
    return std::make_shared<XarrayZarrArray<Store, float>>(dataset.template create_array<float>(
        name, units, "<f4", scale_factor, chunkshape, {"time", "gbxindex"}));
  };

  // create shared pointer to 2-D array in dataset for 0th mass moment in each gridbox over time
  auto m0_ptr = make_array_ptr("massmom0", "", 1.0);

  // const auto c = CombineGDW<Store>{};
  auto m0 = GenericGbxWriter<Store, float, MassMom0Func>(dataset, MassMom0Func{}, m0_ptr, ngbxs);

  return m0;
}

#endif  // LIBS_OBSERVERS2_MASSMOMENTS_WRITERS_HPP_
