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
 * Template for an struct which writes data collected from superdroplets in parallel
 * to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_DO_WRITE_SUPERS_HPP_
#define LIBS_OBSERVERS2_DO_WRITE_SUPERS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "./write_gridbox_to_array.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"

// TODO(CB) attempt 2-D ragged array?

template <typename Store>
struct WriteRaggedCountToArray {
 private:
  std::shared_ptr<XarrayZarrArray<Store, uint32_t>>
      xzarr_ptr;  ///< pointer to time array in dataset

 public:
  WriteRaggedCountToArray(const Dataset<Store> &dataset, const size_t maxchunk)
      : xzarr_ptr(std::make_shared<XarrayZarrArray<Store, uint32_t>>(
            dataset.template create_raggedcount_array<uint32_t>(
                "raggedcount", "", "<u4", 1, {maxchunk}, {"time"}, "superdroplets"))) {}

  /* write_ragged_count_to_array operator() writes the total number of super-droplets in the domain
  "totnsupers" to the raggedcount array in the dataset. Note static conversion from architecture
  dependent, usually 8 byte unsigned integer (size_t = uint64_t), to 4 byte unsigned integer
  (uint32_t). */
  void operator()(const viewd_constsupers totsupers, const Dataset<Store> &dataset) const {
    const auto totnsupers = static_cast<uint32_t>(totsupers.extent(0));
    dataset.write_to_array(xzarr_ptr, totnsupers);
  }

  void write_arrayshape(const Dataset<Store> &dataset) { dataset.write_arrayshape(xzarr_ptr); }
};

/* template class for observer with at_start_step function that collects variables from each
superdroplet in each gridbox in parallel and then writes them to their respective ragged arrays in
a dataset alongside the raggedcount for the arrays */
template <typename Store, WriteGridboxToArray<Store, viewd_constsupers> WriteSupersToArray>
class DoWriteSupers {
 private:
  const Dataset<Store> &dataset;                              ///< dataset to write data to
  WriteRaggedCountToArray<Store> write_raggedcount_to_array;  ///< raggedcount array in dataset
  WriteSupersToArray
      write2array;  ///< object collects superdrops data and writes it to arrays in dataset

  /* Use the writer's functor to collect data from superdroplets in parallel.
  Then write the data to ragged arrays in the dataset */
  void write_superdrops_data(const viewd_constsupers totsupers) const {
    auto functor = write2array.get_functor(totsupers);
    const size_t totnsupers(totsupers.extent(0));
    Kokkos::parallel_for("range_policy_collect_totsupers_data",
                         Kokkos::RangePolicy<ExecSpace>(0, totnsupers), functor);
    write2array.write_to_array(dataset);
  }

  /* Collect data from superdroplets and write into ragged arrays in the dataset alongside the
   * raggedcount */
  void at_start_step(const viewd_constsupers totsupers) const {
    write_superdrops_data(totsupers);
    write_raggedcount_to_array(totsupers, dataset);
  }

 public:
  DoWriteSupers(const Dataset<Store> &dataset, const size_t maxchunk,
                WriteSupersToArray write2array)
      : dataset(dataset),
        write_raggedcount_to_array(WriteRaggedCountToArray(dataset, maxchunk)),
        write2array(write2array) {}

  ~DoWriteSupers() {
    write2array.write_arrayshape(dataset);
    write_raggedcount_to_array.write_arrayshape(dataset);
  }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes write superdrops observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(totsupers);
  }
};

#endif  // LIBS_OBSERVERS2_DO_WRITE_SUPERS_HPP_
