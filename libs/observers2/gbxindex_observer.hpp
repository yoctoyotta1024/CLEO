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
  Dataset<Store> &dataset;                                   ///< dataset to write time data to
  std::shared_ptr<XarrayZarrArray<Store, float>> xzarr_ptr;  ///< pointer to time array in dataset

  // increment size of time dimension in dataset and write out time data to array in the dataset.
  // Note conversion of time from double (8 bytes) to single precision (4 bytes float) in output
  void at_start_step(const unsigned int t_mdl) const {
    const auto ntimes = size_t{dataset.get_dimension("time") + 1};
    const auto timedim = std::pair<std::string, size_t>({"time", ntimes});
    dataset.set_dimension(timedim);

    const auto time = static_cast<float>(step2dimlesstime(t_mdl));
    dataset.write_to_array(xzarr_ptr, time);
  }  // tODO

 public:
  GbxindexObserver(Dataset<Store> &dataset, const size_t maxchunk)
      : dataset(dataset),
        xzarr_ptr(
            std::make_shared<XarrayZarrArray<Store, float>>(dataset.template create_array<float>(
                "time", "s", "<f4", dlc::TIME0, {maxchunk}, {"time"}))) {
    dataset.add_dimension({"time", 0});
  }

  ~GbxindexObserver() { dataset.write_arrayshape(xzarr_ptr); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes gbxindex observer\n";
    // TODO(CB) WIP
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {}

  unsigned int next_obs(const unsigned int t_mdl) const { return LIMITVALUES::uintmax; }

  bool on_step(const unsigned int t_mdl) const { return false; }
};

#endif  // LIBS_OBSERVERS2_GBXINDEX_OBSERVER_HPP_
