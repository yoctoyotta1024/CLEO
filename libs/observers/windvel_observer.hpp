/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: windvel_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 22nd May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to write variables related to Gridboxes' state at the start of
 * a constant interval timestep to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS_WINDVEL_OBSERVER_HPP_
#define LIBS_OBSERVERS_WINDVEL_OBSERVER_HPP_

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
#include "zarr/dataset.hpp"

/**
 * @brief Constructs type sastifying the CollectDataForDataset concept for a given Store (using an
 * instance of the GenericCollectData class) which writes a wind velocity component to an Xarray in
 * a dataset.
 *
 * Function return type writes a wind velocity component to an Xarray as a 4-byte floating point
 * type with units "m/s" by collecting data according to the given FunctorFunc from within a
 * Kokkos::parallel_for loop over gridboxes with a range policy.
 *
 * @param dataset The dataset to write the wind velocity component to.
 * @param ffunc The functor function to collect the wind velocity component from within a parallel
 * range policy over gridboxes.
 * @param maxchunk The maximum chunk size (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return CollectDataForDataset<Store> An instance satisfying the CollectDataForDataset concept for
 * collecting a wind velocity component from each gridbox.
 */
template <typename Store, typename FunctorFunc>
CollectDataForDataset<Store> auto CollectWindVariable(const Dataset<Store> &dataset,
                                                      const FunctorFunc ffunc,
                                                      const std::string_view name,
                                                      const size_t maxchunk, const size_t ngbxs) {
  const auto units = std::string_view("m/s");
  const auto scale_factor = dlc::W0;
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr =
      dataset.template create_array<float>(name, units, scale_factor, chunkshape, dimnames);
  return GenericCollectData(ffunc, xzarr, ngbxs);
}

/**
 * @brief Functor operator to perform a copy of the wvel at the centre of each gridbox "wvel" to
 * d_data within Kokkos::parallel_for loop over gridboxes with range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of wvel from double (8 bytes) to single precision float (4 bytes).
 *
 * @param ii The index of the gridbox.
 * @param d_gbxs The view of gridboxes on device.
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for wvel form each gridbox.
 */
struct WvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto wvel = static_cast<float>(d_gbxs(ii).state.wvelcentre());
    d_data(ii) = wvel;
  }
};

/**
 * @brief Functor operator to perform a copy of the uvel at the centre of each gridbox "uvel" to
 * d_data within Kokkos::parallel_for loop over gridboxes with range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of uvel from double (8 bytes) to single precision float (4 bytes).
 *
 * @param ii The index of the gridbox.
 * @param d_gbxs The view of gridboxes on device.
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for uvel form each gridbox.
 */
struct UvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto uvel = static_cast<float>(d_gbxs(ii).state.uvelcentre());
    d_data(ii) = uvel;
  }
};

/**
 * @brief Functor operator to perform a copy of the vvel at the centre of each gridbox "vvel" to
 * d_data within Kokkos::parallel_for loop over gridboxes with range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of vvel from double (8 bytes) to single precision float (4 bytes).
 *
 * @param ii The index of the gridbox.
 * @param d_gbxs The view of gridboxes on device.
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for vvel form each gridbox.
 */
struct VvelFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto vvel = static_cast<float>(d_gbxs(ii).state.vvelcentre());
    d_data(ii) = vvel;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept for collecting all three
 * wind velocity components in each gridbox and writing them to a dataset.
 *
 * This function combines CollectDataForDataset types for the three wind velocity components (wvel,
 * vvel and uvel) in each gridbox using instances of the GenericCollectData class.
 *
 * @tparam Store The type of the dataset store.
 * @param dataset The dataset to write the wind velocity components to.
 * @param maxchunk The maximum chunk size (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return CollectDataForDataset<Store> An instance of CollectDataForDataset for collecting wind
 * velocity data.
 */
template <typename Store>
inline CollectDataForDataset<Store> auto CollectWindVel(const Dataset<Store> &dataset,
                                                        const size_t maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Store> auto wvel =
      CollectWindVariable<Store, WvelFunc>(dataset, WvelFunc{}, "wvel", maxchunk, ngbxs);

  const CollectDataForDataset<Store> auto uvel =
      CollectWindVariable<Store, UvelFunc>(dataset, UvelFunc{}, "uvel", maxchunk, ngbxs);

  const CollectDataForDataset<Store> auto vvel =
      CollectWindVariable<Store, VvelFunc>(dataset, VvelFunc{}, "vvel", maxchunk, ngbxs);

  return vvel >> uvel >> wvel;
}

/**
 * @brief Constructs an observer which writes the wind velocity components in each gridbox (wvel,
 * vvel and uvel) at start of each observation timestep to a arrays with a constant observation
 * timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Observer An observer instance for writing the wind velocity components.
 */
template <typename Store>
inline Observer auto WindVelObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                     const size_t maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Store> auto windvel = CollectWindVel(dataset, maxchunk, ngbxs);
  return WriteToDatasetObserver(interval, dataset, windvel);
}

#endif  // LIBS_OBSERVERS_WINDVEL_OBSERVER_HPP_
