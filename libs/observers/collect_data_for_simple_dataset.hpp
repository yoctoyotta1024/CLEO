/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: collect_data_for_simple_dataset.hpp
 * Project: observers
 * Created Date: Tuesday 15th July 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Operator to combine types which satisfy the CollectDataForDataset concept when the dataset
 * is a SimpleDataset with a FSStore. Operator useful e.g. to make construction of various observers
 * easier when combining multiple "CollectData" types which satisfy the
 * CollectDataForDataset<SimpleDataset<FSStore>> concept.
 */

#ifndef LIBS_OBSERVERS_COLLECT_DATA_FOR_SIMPLE_DATASET_HPP_
#define LIBS_OBSERVERS_COLLECT_DATA_FOR_SIMPLE_DATASET_HPP_

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "observers/collect_data_for_dataset.hpp"
#include "zarr/fsstore.hpp"
#include "zarr/simple_dataset.hpp"

/**
 * @brief Overloaded operator >> to combine two CollectDataForDataset instances into a new one.
 *
 * @param a First CollectDataForDataset with Dataset=SimpleDataset<FSStore>.
 * @param b Second CollectDataForDataset with Dataset=SimpleDataset<FSStore>.
 * @return CombinedCollectDataForDataset<Obs1, Obs2> Combined CollectDataForDataset.
 */
template <CollectDataForDataset<SimpleDataset<FSStore>> CollectData1,
          CollectDataForDataset<SimpleDataset<FSStore>> CollectData2>
auto operator>>(const CollectData1 a, const CollectData2 b) {
  return CombinedCollectDataForDataset<CollectData1, CollectData2>(a, b);
}

#endif  // LIBS_OBSERVERS_COLLECT_DATA_FOR_SIMPLE_DATASET_HPP_
