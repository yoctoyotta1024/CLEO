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
 * Last Modified: Friday 29th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Template for an struct which writes data collected from Gridboxes in parallel
 * to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_
#define LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "./write_gridboxes.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"

/* template class for observer with at_start_step function that collects variables from each
gridbox in parallel and then writes them to their repspective arrays in a dataset */
template <typename Store, GridboxDataWriter<Store> GbxWriter>
class WriteGridboxes {
 private:
  Dataset<Store> &dataset;  ///< dataset to write data to
  GbxWriter writer;  ///< object collects data from gridboxes and writes it to arrays in the dataset

  // use functor from writer to collect data from gridboxes in parallel
  void collect_data_from_gridboxes(const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    auto functor = writer.get_functor(d_gbxs);
    Kokkos::parallel_for("collect_gbxs_data", Kokkos::RangePolicy<ExecSpace>(0, ngbxs), functor);
  }

  // collect data from gridboxes and then write it to arrays in the dataset
  void at_start_step(const viewd_constgbx d_gbxs) const {
    collect_data_from_gridboxes(d_gbxs);
    writer.write_to_array(dataset);
  }

 public:
  WriteGridboxes(Dataset<Store> &dataset, GbxWriter writer) : dataset(dataset), writer(writer) {}

  ~WriteGridboxes() { writer.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes write gridboxes observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(d_gbxs);
  }
};

#endif  // LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_
