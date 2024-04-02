/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: do_write_supers.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 2nd April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Template for an struct which writes data collected from super-droplets in parallel
 * to individual ragged arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_DO_WRITE_SUPERS_HPP_
#define LIBS_OBSERVERS2_DO_WRITE_SUPERS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "./nsupers_observer.hpp"
#include "./write_gridbox_to_array.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"

/* template class for observer with at_start_step function that collects variables from all the
super-droplets in each gridbox in parallel and then writes them to their respective ragged arrays in
a dataset */
template <typename Store, typename ParallelLoopPolicy,
          WriteGridboxToArray<Store> WriteSupersToArray>
class DoWriteSupers {
 private:
  const Dataset<Store> &dataset;     ///< dataset to write data to
  WriteSupersToArray write2array;    ///< object collects superdrops data and writes it to arrays
                                     ///< in dataset, including array for raggedcount
  ParallelLoopPolicy parallel_loop;  ///< function like object to call during at_start_step to
                                     ///< loop over gridboxes

  /* Use the writer's functor to collect data from gridboxes in parallel.
  Then write the datat to arrays in the dataset */
  void at_start_step(const viewd_constgbx d_gbxs) const {
    auto functor = write2array.get_functor(d_gbxs);
    parallel_loop(functor, d_gbxs);
    write2array.write_to_array(dataset);
  }

 public:
  DoWriteSupers(ParallelLoopPolicy parallel_loop, const Dataset<Store> &dataset,
                WriteSupersToArray write2array)
      : dataset(dataset), write2array(write2array), parallel_loop(parallel_loop) {}

  ~DoWriteSupers() { write2array.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes write supers observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs) const {
    at_start_step(d_gbxs);
  }
};

#endif  // LIBS_OBSERVERS2_DO_WRITE_SUPERS_HPP_
