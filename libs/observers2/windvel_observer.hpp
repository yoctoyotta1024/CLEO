/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: windvel_observer.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 4th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to write variables related to Gridboxes' state at the start of
 * a constant interval timestep to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_WINDVEL_OBSERVER_HPP_
#define LIBS_OBSERVERS2_WINDVEL_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./collect_data_for_dataset.hpp"
#include "./generic_collect_data.hpp"
#include "./write_to_dataset_observer.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"

/* returns CollectDataForDataset which writes a state variable from
each gridbox to an array in a dataset in a given store for a given datatype and using a given
function-like functor */
template <typename Store, typename FunctorFunc>
CollectDataForDataset<Store> auto CollectWindVariable(const Dataset<Store> &dataset,
                                                      const FunctorFunc ffunc,
                                                      const std::string_view name,
                                                      const size_t maxchunk, const size_t ngbxs) {
  const auto units = std::string_view("m/s");
  const auto dtype = std::string_view("<f4");
  const auto scale_factor = dlc::W0;
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr =
      dataset.template create_array<float>(name, units, dtype, scale_factor, chunkshape, dimnames);
  return GenericCollectData(ffunc, xzarr, ngbxs);
}

/* Operator is functor to perform copy of wvel at the centre of each gridbox to d_data
in parallel. Note conversion of wvel from double (8 bytes) to single precision (4 bytes
float) in output */
struct WvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto wvel = static_cast<float>(d_gbxs(ii).state.wvelcentre());
    d_data(ii) = wvel;
  }
};

/* Operator is functor to perform copy of uvel at the centre of each gridbox to d_data
in parallel. Note conversion of uvel from double (8 bytes) to single precision (4 bytes
float) in output */
struct UvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto uvel = static_cast<float>(d_gbxs(ii).state.uvelcentre());
    d_data(ii) = uvel;
  }
};

/* Operator is functor to perform copy of vvel at the centre of each gridbox to d_data
in parallel. Note conversion of vvel from double (8 bytes) to single precision (4 bytes
float) in output */
struct VvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto vvel = static_cast<float>(d_gbxs(ii).state.vvelcentre());
    d_data(ii) = vvel;
  }
};

/* constructs CollectDataForDataset for a given Store which writes the wind velocity at center of
each gridbox using an instance of the GenericCollectData class */
template <typename Store>
inline CollectDataForDataset<Store> auto CollectWindVel(const Dataset<Store> &dataset,
                                                        const int maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Store> auto wvel =
      CollectWindVariable<Store, WvelFunc>(dataset, WvelFunc{}, "wvel", maxchunk, ngbxs);

  const CollectDataForDataset<Store> auto uvel =
      CollectWindVariable<Store, UvelFunc>(dataset, UvelFunc{}, "uvel", maxchunk, ngbxs);

  const CollectDataForDataset<Store> auto vvel =
      CollectWindVariable<Store, VvelFunc>(dataset, VvelFunc{}, "vvel", maxchunk, ngbxs);

  return vvel >> uvel >> wvel;
}

/* constructs observer which writes writes the wind velocity at center of
each gridbox with a constant timestep 'interval' using an instance of the WriteToDatasetObserver
class */
template <typename Store>
inline Observer auto WindVelObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                     const int maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Store> auto windvel = CollectWindVel(dataset, maxchunk, ngbxs);
  return WriteToDatasetObserver(interval, dataset, windvel);
}

#endif  // LIBS_OBSERVERS2_WINDVEL_OBSERVER_HPP_
