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

/* Operator is functor to perform copy of number of superdrops in each gridbox to d_data
  in parallel. Note conversion of nsupers from size_t (8 bytes) to single precision (4 bytes
  unsigned integer) in output */
struct RaggedCountFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto nsupers = static_cast<uint32_t>(d_gbxs(ii).supersingbx.nsupers());
    d_data(ii) = nsupers;
  }
};

/* returns WriteGridboxToArray which writes the number of superdrops in each
gridbox to an array in a dataset in a store called "raggedcount_nsupers" */
template <typename Store>
WriteGridboxToArray<Store> auto RaggedCountWriter(const Dataset<Store> &dataset,
                                                  const size_t maxchunk, const size_t ngbxs) {
  return GenericWriteSupersToXarray<Store, uint32_t, RaggedCountFunc, XarrayForSupersRaggedCount>(
      dataset, "raggedcount_nsupers", "<u4", maxchunk, ngbxs, RaggedCountFunc{});
}

/* template class for observer with at_start_step function that collects variables from each
superdroplet in each gridbox in parallel and then writes them to their repspective ragged arrays in
a dataset alongside the raggedcount for the arrays */
template <typename Store, WriteGridboxToArray<Store> WriteSupersToArray,
          WriteGridboxToArray<Store> WriteRaggedCountToArray>
class DoWriteSupers {
 private:
  const Dataset<Store> &dataset;  ///< dataset to write data to
  WriteSupersToArray
      write2array;  ///< object collects superdrops data and writes it to arrays in dataset
  WriteRaggedCountToArray write2raggedcount;  ///< raggedcount array in dataset

  /* Use the writer's functor to collect data from superdroplets in parallel.
  Then write the data to ragged arrays in the dataset */
  void write_superdrops_data(const viewd_constsupers totsupers) const {
    auto functor = write2array.get_functor(totsupers);
    const size_t totnsupers(totsupers.extent(0));
    Kokkos::parallel_for("range_policy_collect_totsupers_data",
                         Kokkos::RangePolicy<ExecSpace>(0, totnsupers), functor);
    write2array.write_to_array(dataset);
  }

  /* Use the raggedcount writer's functor to collect the raggedcount (i.e. the number of
  superdroplets in each gridbox) and write it to an array in the dataset */
  void write_raggedcount(const viewd_constgbx d_gbxs) const {
    auto functor = write2raggedcount.get_functor(d_gbxs);
    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for("range_policy_collect_gbxs_data", Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
                         functor);
    write2raggedcount.write_to_array(dataset);
  }

  /* Collect data from superdroplets and write into ragged arrays in the dataset alongside the
   * raggedcount */
  void at_start_step(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) const {
    write_superdrops_data(totsupers);
    write_raggedcount(d_gbxs);
  }

 public:
  DoWriteSupers(const Dataset<Store> &dataset, WriteSupersToArray write2array)
      : dataset(dataset),
        write2array(write2array),
        write2raggedcount(RaggedCountWriter(dataset, maxchunk, ngbxs)) {}

  ~DoWriteSupers() { write2array.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes write superdrops observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(d_gbxs, totsupers);
  }
};

#endif  // LIBS_OBSERVERS2_DO_WRITE_SUPERS_HPP_
