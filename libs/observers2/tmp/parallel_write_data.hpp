/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: parallel_write_data.hpp
 * Project: tmp
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 3rd April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Template for "ParallelWriteData" function-like object (see write_to_dataset_observer.hpp) for
 * writing data from gridboxes and/or superdroplets to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_TMP_PARALLEL_WRITE_DATA_HPP_
#define LIBS_OBSERVERS2_TMP_PARALLEL_WRITE_DATA_HPP_

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "zarr2/dataset.hpp"

/* struct for "ParallelWriteData" (see write_to_dataset_observer.hpp) to collect data from
gridboxes in a parallel loop and write it to arrays in a dataset */
template <typename Store, typename CollectData>
class ParallelWriteGridboxes {
 private:
  const Dataset<Store> &dataset;  ///< dataset to write data to
  CollectData collect_data;  ///< functions to collect data within gbxs loop and write in dataset

 public:
  ParallelWriteGridboxes(const Dataset<Store> &dataset, CollectData collect_data)
      : dataset(dataset), collect_data(collect_data) {}

  ~ParallelWriteGridboxes() { collect_data.write_arrayshapes(dataset); }

  /* Use the writer's functor to collect data from gridboxes in parallel.
    Then write the data in the dataset */
  void operator()(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) {
    auto functor = collect_data.get_functor(d_gbxs);

    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for("write_gridboxes", Kokkos::RangePolicy<ExecSpace>(0, ngbxs), functor);

    collect_data.write_to_arrays(dataset);
  }
};

/* struct for "ParallelWriteData" (see write_to_dataset_observer.hpp) to collect data from
superdroplets in a parallel loop and write it to ragged arrays in a dataset */
template <typename Store, typename CollectData, typename RaggedCount>
class ParallelWriteSupers {
 private:
  const Dataset<Store> &dataset;  ///< dataset to write data to
  CollectData collect_data;  ///< functions to collect data within supers loop and write in dataset
  RaggedCount ragged_count;  ///< functions to write ragged count variable in dataset

 public:
  ParallelWriteSupers(const Dataset<Store> &dataset, CollectData collect_data,
                      RaggedCount write_ragged_count)
      : dataset(dataset), collect_data(collect_data), write_ragged_count(write_ragged_count) {}

  ~ParallelWriteSupers() {
    collect_data.write_arrayshapes(dataset);
    ragged_count.write_arrayshape(dataset);
  }

  /* Use the writer's functor to collect data from superdrops in parallel.
  Then write the data in the dataset */
  void operator()(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) {
    auto functor = collect_data.get_functor(totsupers);

    const size_t totnsupers(totsupers.extent(0));
    Kokkos::parallel_for("write_supers", Kokkos::RangePolicy<ExecSpace>(0, totnsupers), functor);

    collect_data.write_to_arrays(dataset);

    ragged_count.write_to_array(dataset, totsupers);
  }
};

#endif  // LIBS_OBSERVERS2_TMP_PARALLEL_WRITE_DATA_HPP_
