/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: collect_data_for_collective_dataset.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 15th July 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Operator to combine types which satisfy the CollectDataForDataset concept when the dataset
 * is a CollectiveDataset with a FSStore and CartesianDecomposition. Operator useful e.g. to make
 * construction of various observers easier when combining multiple "CollectData" types which
 * satisfy the CollectDataForDataset<CollectiveDataset<FSStore, CartesianDecomposition>> concept.
 */

#ifndef LIBS_CARTESIANDOMAIN_COLLECT_DATA_FOR_COLLECTIVE_DATASET_HPP_
#define LIBS_CARTESIANDOMAIN_COLLECT_DATA_FOR_COLLECTIVE_DATASET_HPP_

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesian_decomposition.hpp"
#include "observers/collect_data_for_dataset.hpp"
#include "zarr/collective_dataset.hpp"
#include "zarr/fsstore.hpp"

/**
 * @brief Overloaded operator >> to combine two CollectDataForDataset instances into a new one.
 *
 * @param a First CollectDataForDataset with
 *            Dataset=CollectiveDataset<FSStore, CartesianDecomposition>.
 * @param b Second CollectDataForDataset with
 *            Dataset=CollectiveDataset<FSStore, CartesianDecomposition>.
 * @return CombinedCollectDataForDataset<Obs1, Obs2> Combined CollectDataForDataset.
 */
template <CollectDataForDataset<CollectiveDataset<FSStore, CartesianDecomposition>> CollectData1,
          CollectDataForDataset<CollectiveDataset<FSStore, CartesianDecomposition>> CollectData2>
auto operator>>(const CollectData1 a, const CollectData2 b) {
  return CombinedCollectDataForDataset<CollectData1, CollectData2>(a, b);
}

#endif  // LIBS_CARTESIANDOMAIN_COLLECT_DATA_FOR_COLLECTIVE_DATASET_HPP_
