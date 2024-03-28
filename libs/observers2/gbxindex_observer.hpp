/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: gbxindex_observer.hpp
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
 * Observer to output gridbox indexes at the start of each simulation to an array in a dataset
 */

#ifndef LIBS_OBSERVERS2_GBXINDEX_OBSERVER_HPP_
#define LIBS_OBSERVERS2_GBXINDEX_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <functional>
#include <iostream>
#include <memory>
#include <utility>

#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"

/* template observer which writes gbxindex for every gridbox out to a 1-D array
as a coordinate of an xarray dataset */
template <typename Store>
class GbxindexObserver {
 private:
  Dataset<Store> &dataset;                                    ///< dataset to write gbxindex data to
  std::shared_ptr<XarrayZarrArray<Store, size_t>> xzarr_ptr;  ///< pointer to gbxindex array

  /* increment size of gbxindex dimension in dataset and write out gbxindex data
  to array in the dataset. */
  void write_gbxindex(const viewd_constgbx d_gbxs) const {
    // dataset.write_to_array(xzarr_ptr, h_data);
  }

 public:
  GbxindexObserver(Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs)
      : dataset(dataset),
        xzarr_ptr(std::make_shared<XarrayZarrArray<Store, size_t>>(
            dataset.template create_coordinate_array<size_t>("gbxindex", "", "<u4", 1, maxchunk,
                                                             ngbxs))) {}

  ~GbxindexObserver() { dataset.write_arrayshape(xzarr_ptr); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes gbxindex observer\n";
    write_gbxindex(d_gbxs);
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {}

  unsigned int next_obs(const unsigned int t_mdl) const { return LIMITVALUES::uintmax; }

  bool on_step(const unsigned int t_mdl) const { return false; }
};

#endif  // LIBS_OBSERVERS2_GBXINDEX_OBSERVER_HPP_
