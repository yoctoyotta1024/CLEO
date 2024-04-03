/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: parallel_write_data.hpp
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
 * Template for "ParallelWriteData" function-like object (see write_to_dataset_observer.hpp) for
 * writing data from gridboxes and/or superdroplets to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_PARALLEL_WRITE_DATA_HPP_
#define LIBS_OBSERVERS2_PARALLEL_WRITE_DATA_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../kokkosaliases.hpp"
#include "zarr2/dataset.hpp"

/* struct for "ParallelWriteData" (see write_to_dataset_observer.hpp) to collect data from
gridboxes in a parallel loop and write it to arrays in a dataset */
template <typename Store, CollectDataForDataset<Store> CollectData>
class ParallelWriteGridboxes {
 private:
  const Dataset<Store> &dataset;  ///< dataset to write data to
  CollectData collect_data;  ///< functions to collect data within gbxs loop and write in dataset

  /* parallel loop over gridboxes using Kokkos Range Policy */
  void parallel_write_gridboxes(const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for("write_gridboxes", Kokkos::RangePolicy<ExecSpace>(0, ngbxs), functor);
  }

 public:
  ParallelWriteGridboxes(const Dataset<Store> &dataset, CollectData collect_data)
      : dataset(dataset), collect_data(collect_data) {}

  ~ParallelWriteGridboxes() { collect_data.write_arrayshapes(dataset); }

  /* Use the CollectData instance's functor to collect data from gridboxes in a parallel loop.
    Then write the data in the dataset. Inclusion of totsupers so that object can be used as
    "ParallelWriteData" function in DoWriteToDataset struct */
  void operator()(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) const {
    auto functor = collect_data.get_functor(d_gbxs, totsupers);
    parallel_write_gridboxes(d_gbxs);
    collect_data.write_to_arrays(dataset);
  }
};

/**
 * @brief Concept for CollectRaggedCount is all types that have functions for writing the ragged
 * count of superdroplet arrays to an array in a dataset.
 *
 * @tparam CRC The type that satisfies the CollectRaggedCount concept.
 */
template <typename CRC, typename Store>
concept CollectRaggedCount =
    requires(CRC crc, const Dataset<Store> &ds, const viewd_constsupers totsupers) {
      { crc.write_to_array(ds, totsupers) } -> std::same_as<void>;
      { crc.write_arrayshape(ds) } -> std::same_as<void>;
    };

/* struct for "ParallelWriteData" (see write_to_dataset_observer.hpp) to collect data from
superdroplets in a parallel loop and write it to ragged arrays in a dataset */
template <typename Store, CollectDataForDataset<Store> CollectData,
          CollectRaggedCount<Store> RaggedCount>
class ParallelWriteSupers {
 private:
  const Dataset<Store> &dataset;  ///< dataset to write data to
  CollectData collect_data;  ///< functions to collect data within supers loop and write in dataset
  RaggedCount ragged_count;  ///< functions to write ragged count variable in dataset

  /* parallel loop over superdroplets using Kokkos Range Policy */
  void parallel_write_supers(const viewd_constsupers totsupers) const {
    const size_t totnsupers(totsupers.extent(0));
    Kokkos::parallel_for("write_supers", Kokkos::RangePolicy<ExecSpace>(0, totnsupers), functor);
  }

 public:
  ParallelWriteSupers(const Dataset<Store> &dataset, CollectData collect_data,
                      RaggedCount ragged_count)
      : dataset(dataset), collect_data(collect_data), ragged_count(ragged_count) {}

  ~ParallelWriteSupers() {
    collect_data.write_arrayshapes(dataset);
    ragged_count.write_arrayshape(dataset);
  }

  /* Use the CollectData instance's functor to collect data from superdroplets in a parallel loop.
    Then write the data in the dataset alongside ragged count for arrays. Inclusion of d_gbxs
    so that object can be used as "ParallelWriteData" function in DoWriteToDataset struct */
  void operator()(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) const {
    collect_data.reallocate_views(totsupers.extent(0));
    auto functor = collect_data.get_functor(d_gbxs, totsupers);
    parallel_write_supers(totsupers);
    collect_data.write_to_arrays(dataset);
    ragged_count.write_to_array(dataset, totsupers);
  }
};

#endif  // LIBS_OBSERVERS2_PARALLEL_WRITE_DATA_HPP_
