/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: totnsupers_observer.hpp
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
 * Observer to output total number of superdroplets at the start of each timestep to an array in a
 * dataset
 */

#ifndef LIBS_OBSERVERS2_TOTNSUPERS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_TOTNSUPERS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/buffer.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"

/* template observer which writes gbxindex for every gridbox out to a 1-D array
as a coordinate of an xarray dataset */
template <typename Store>
class DoTotNsupersObs {
 private:
  Dataset<Store> &dataset;  ///< dataset to write totnsupers data to
  std::shared_ptr<XarrayZarrArray<Store, uint32_t>> xzarr_ptr;  ///< pointer to totnsupers array

  void at_start_step(const viewd_constsupers totsupers) const {
    const auto data = static_cast<uint32_t>(totsupers.extent(0));
    dataset.write_to_array(xzarr_ptr, data);
  }

 public:
  DoTotNsupersObs(Dataset<Store> &dataset, const size_t maxchunk)
      : dataset(dataset),
        xzarr_ptr(std::make_shared<XarrayZarrArray<Store, uint32_t>>(
            dataset.template create_array<uint32_t>("totnsupers", "", "<u4", 1, {maxchunk},
                                                    {"time"}))) {}

  ~DoTotNsupersObs() { dataset.write_arrayshape(xzarr_ptr); }

  /* write the gbxindex of every gridbox in d_gbxs to the gbxindex array in the dataset and assert
  the size of the gbxindex dimension in the dataset is correct */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes totnsupers observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(totsupers);
  }
};

/* constructs observer which writes time variable out to a 1-D array with a constant timestep
'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto TotNsupersObserver(const unsigned int interval, Dataset<Store> &dataset,
                                        const int maxchunk) {
  return ConstTstepObserver(interval, DoTotNsupersObs(dataset, maxchunk));
}

#endif  // LIBS_OBSERVERS2_TOTNSUPERS_OBSERVER_HPP_
