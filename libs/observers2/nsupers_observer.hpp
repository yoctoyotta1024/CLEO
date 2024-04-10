/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: nsupers_observer.hpp
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
 * Observer to write the number of superdroplets in each gridbox at the start of
 * a constant interval timestep to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_NSUPERS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_NSUPERS_OBSERVER_HPP_

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

/* Operator is functor to perform copy of the number of superdroplets in each gridbox "nsupers" to
d_data in parallel. Note conversion of nsupers from size_t (architecture dependent usually long
unsigned int = 8 bytes) to single precision (uint32_t = 4 bytes) in output */
struct NsupersFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto nsupers = static_cast<uint32_t>(d_gbxs(ii).supersingbx.nsupers());
    d_data(ii) = nsupers;
  }
};

/* constructs CollectDataForDataset for a given Store which writes the number of superdroplets
in each gridbox using an instance of the GenericCollectData class */
template <typename Store>
inline CollectDataForDataset<Store> auto CollectNsupers(const Dataset<Store> &dataset,
                                                        const int maxchunk, const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr =
      dataset.template create_array<uint32_t>("nsupers", "", "<u4", 1, chunkshape, dimnames);
  return GenericCollectData(NsupersFunc{}, xzarr, ngbxs);
}

/* constructs observer which writes the number of superdroplets each gridbox
with a constant timestep 'interval' using an instance of the WriteToDatasetObserver class */
template <typename Store>
inline Observer auto NsupersObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                     const int maxchunk, const size_t ngbxs) {
  return WriteToDatasetObserver(interval, dataset, CollectNsupers(dataset, maxchunk, ngbxs));
}

#endif  // LIBS_OBSERVERS2_NSUPERS_OBSERVER_HPP_
