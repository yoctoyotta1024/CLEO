/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: nsupers_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to write the number of superdroplets in each gridbox at the start of
 * a constant interval timestep to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS_NSUPERS_OBSERVER_HPP_
#define LIBS_OBSERVERS_NSUPERS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "observers/collect_data_for_dataset.hpp"
#include "observers/generic_collect_data.hpp"
#include "observers/write_to_dataset_observer.hpp"

/**
 * @brief Functor operator to perform a copy of the number of superdroplets in each gridbox
 * "nsupers" to d_data within Kokkos::parallel_for loop over gridboxes with range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of nsupers from size_t (architecture-dependent, usually long unsigned int = 8
 * bytes) to single precision (uint32_t = 4 bytes).
 *
 * @param ii The index of the gridbox.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the number of superdroplets.
 */
struct NsupersFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto nsupers = static_cast<uint32_t>(d_gbxs(ii).supersingbx.nsupers());
    d_data(ii) = nsupers;
  }
};

/**
 * @brief Constructs type sastifying the CollectDataForDataset concept for a given Dataset (using an
 * instance of the GenericCollectData class) which writes the number of superdroplets in each
 * gridbox "nsupers" during the functor call.
 *
 * @tparam Dataset Type of dataset
 * @param dataset The dataset to write nsupers to.
 * @param maxchunk The maximum chunk size (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return CollectDataForDataset<Dataset> An instance satisfying the CollectDataForDataset concept
 * for collecting the number of superdroplets in each gridbox.
 */
template <typename Dataset>
inline CollectDataForDataset<Dataset> auto CollectNsupers(const Dataset &dataset,
                                                          const size_t maxchunk,
                                                          const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr =
      dataset.template create_array<uint32_t>("nsupers", "", 1, chunkshape, dimnames);
  return GenericCollectData(NsupersFunc{}, xzarr, ngbxs);
}

/**
 * @brief Constructs an observer which writes the number of superdroplets in each gridbox "nsupers"
 * at start of each observation timestep to an array with a constant observation timestep
 * "interval".
 *
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Observer An observer instance for writing the number of superdroplets.
 */
template <typename Dataset>
inline Observer auto NsupersObserver(const unsigned int interval, const Dataset &dataset,
                                     const size_t maxchunk, const size_t ngbxs) {
  return WriteToDatasetObserver(interval, dataset, CollectNsupers(dataset, maxchunk, ngbxs));
}

#endif  // LIBS_OBSERVERS_NSUPERS_OBSERVER_HPP_
