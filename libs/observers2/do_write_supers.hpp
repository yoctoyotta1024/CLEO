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
 * Last Modified: Wednesday 3rd April 2024
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

struct ParallelSupersRangePolicy {
  template <typename Functor>
  void operator()(Functor functor, const viewd_constsupers totsupers) const {
    const size_t totnsupers(totsupers.extent(0));
    Kokkos::parallel_for("range_policy_collect_totsupers_data",
                         Kokkos::RangePolicy<ExecSpace>(0, ntotsupers), functor);
  }
};

/* template class for observer with at_start_step function that collects variables from all the
super-droplets in each gridbox in parallel and then writes them to their respective ragged arrays in
a dataset */
template <typename Store, WriteGridboxToArray<Store> WriteSupersToArray>
class DoWriteSupers {
 private:
  ParallelSupersRangePolicy parallel_loop;  ///< function like object to call during at_start_step
                                            ///< to loop over gridboxes
  const Dataset<Store> &dataset;            ///< dataset to write data to
  WriteSupersToArray write2array;  ///< object collects superdrops data and writes it to arrays
                                   ///< in dataset, including array for raggedcount

  /* Use the writer's functor to collect data from superdroplets in parallel.
  Then write the data to arrays in the dataset */
  void at_start_step(const viewd_constsupers totsupers) const {
    auto functor = write2array.get_functor(totsupers);
    parallel_loop(functor, totsupers);
    write2array.write_to_array(dataset);
  }

 public:
  DoWriteSupers(const Dataset<Store> &dataset, WriteSupersToArray write2array)
      : parallel_loop(ParallelSupersRangePolicy{}), dataset(dataset), write2array(write2array) {}

  ~DoWriteSupers() { write2array.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes write supers observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(totsupers);
  }
};

#endif  // LIBS_OBSERVERS2_DO_WRITE_SUPERS_HPP_
