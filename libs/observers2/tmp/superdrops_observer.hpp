/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: superdrops_observer.hpp
 * Project: tmp
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
 * Observer to write variables related to superdroplet attributes at the start of
 * a constant interval timestep to ragged arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_TMP_SUPERDROPS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_TMP_SUPERDROPS_OBSERVER_HPP_

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
#include "superdrops/superdrop.hpp"
#include "zarr2/dataset.hpp"

/* returns CollectDataForDataset which writes a variable (e.g. an attribute) from
each superdroplet to an array in a dataset in a given store for a given datatype and using a given
function-like functor */
template <typename Store, typename T, typename FunctorFunc>
CollectDataForDataset<Store> auto CollectSuperdropVariable(
    const Functorfunc ffunc, const std::string_view name, const std::string_view units,
    const std::string_view dtype, const double scale_factor, const size_t maxchunk) {
  const auto chunkshape = std::vector<size_t>{maxchunk};
  const auto dimnames = std::vector<std::string>{"time"};
  const auto sampledimname = std::string_view("superdroplets");
  const auto xzarr_ptr = std::make_shared(dataset.template create_raggedarray<T>(
      name, units, dtype, scale_factor, chunkshape, dimnames, sampledimname));

  return GenericCollectData(ffunc, xzarr_ptr, 0);
}

/* constructs observer which writes writes superdroplet variables (e.g. an attributes) from each
superdroplet with a constant timestep 'interval' using an instance of the WriteToDatasetObserver
class */
template <typename Store>
inline Observer auto SuperdropsObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                        const int maxchunk) {
  const CollectDataForDataset<Store> auto xi =
      CollectSuperdropVariable<Store, uint32_t, XiFunc>(XiFunc{}, "xi", "", "<u8", 1, maxchunk);

  const CollectRaggedCount<Store> auto ragged_count = 1;

  return WriteToDatasetObserver(interval, dataset, collect_superdropsdata, ragged_count);
}

#endif  // LIBS_OBSERVERS2_TMP_SUPERDROPS_OBSERVER_HPP_
