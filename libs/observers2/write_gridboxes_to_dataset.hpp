/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: write_gridboxes_to_dataset.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 28th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Template for an observer which outputs variables from Gridboxes at the start of
 * each timestep in parallel to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_WRITE_GRIDBOXES_TO_DATASET_HPP_
#define LIBS_OBSERVERS2_WRITE_GRIDBOXES_TO_DATASET_HPP_

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"

/* template for observing variables from each gridbox in parallel and then writing them to
their repspective arrays in a dataset */
template <typename Store, typename GridboxesToArrays>
class WriteGridboxesToDataset {
 private:
  Dataset<Store> &dataset;
  GridboxesToArrays gbxs2arrays;

  void collect_data_from_gridboxes(const viewd_constgbx d_gbxs) const {
    auto functor = gbxs2arrays.get_functor(d_gbxs);

    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for("gbxs2arrays", Kokkos::RangePolicy<ExecSpace>(0, ngbxs), functor);
  }

  void at_start_step(const viewd_constgbx d_gbxs) const {
    collect_data_from_gridboxes(d_gbxs);
    gbxs2arrays.write_data(dataset);

    // dataset.set_dimension({"time", time+1}); // TODO(CB) do this with coord observer
  }

 public:
  WriteGridboxesToDataset(Dataset<Store> &dataset, GridboxesToArrays gbxs2arrays)
      : dataset(dataset), gbxs2arrays(gbxs2arrays) {}

  ~WriteGridboxesToDataset() { gbxs2arrays.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes Gridboxes to Dataset observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(d_gbxs);
  }
};

#endif  // LIBS_OBSERVERS2_WRITE_GRIDBOXES_TO_DATASET_HPP_
