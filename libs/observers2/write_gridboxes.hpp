/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: write_gridboxes.hpp
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
 * Template for an observer which writes data collected from Gridboxes at the start of
 * each timestep in parallel to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_
#define LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"

/**
 * @brief Concept for GridboxDataWriter is all types that have functions for creating a functor to
 * collect data from gridboxes and then writing that data to an array in a dataset.
 *
 * @tparam GDW The type that satisfies the GridboxDataWriter concept.
 */
template <typename GDW, typename Store>
concept GridboxDataWriter = requires(GDW gdw, Dataset<Store> &ds, const viewd_constgbx d_gbxs) {
  { gdw.get_functor(d_gbxs) };
  { gdw.write_to_array(ds) } -> std::same_as<void>;
  { gdw.write_arrayshape(ds) } -> std::same_as<void>;
};

/* template class for observing variables from each gridbox in parallel
and then writing them to their repspective arrays in a dataset */
template <typename Store, GridboxDataWriter<Store> GbxWriter>
class WriteGridboxes {
 private:
  Dataset<Store> &dataset;  ///< dataset to write data to
  GbxWriter writer;  ///< object that collects data ffrom girdboxes and writes it to the dataset

  void collect_data_from_gridboxes(const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    auto functor = writer.get_functor(d_gbxs);
    Kokkos::parallel_for("collect_gbxs_data", Kokkos::RangePolicy<ExecSpace>(0, ngbxs), functor);
  }

  void at_start_step(const viewd_constgbx d_gbxs) const {
    collect_data_from_gridboxes(d_gbxs);
    writer.write_to_array(dataset);
    // dataset.set_dimension({"time", time+1}); // TODO(CB) do this with coord observer
  }

 public:
  WriteGridboxes(Dataset<Store> &dataset, GbxWriter writer) : dataset(dataset), writer(writer) {}

  ~WriteGridboxes() { writer.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes Gridboxes to Dataset observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(d_gbxs);
  }
};

#endif  // LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_
