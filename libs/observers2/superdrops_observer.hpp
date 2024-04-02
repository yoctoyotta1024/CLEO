/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: superdrops_observer.hpp
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
 * Observer to output variables related to Gridboxes' state at the start of
 * each timestep to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <array>
#include <concepts>
#include <iostream>
#include <memory>
#include <string_view>

#include "../cleoconstants.hpp"
#include "./observers.hpp"
#include "./write_gridbox_to_array.hpp"
#include "zarr2/buffer.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"
#include "zarr2/zarr_array.hpp"

/* template class for observing superdroplets' attributes in each gridbox to ragged arrays in a
dataset in a store */
template <typename Store, WriteGridboxToArray<Store> WriteGbx>
class DoSuperdropsObs {
 private:
  Dataset<Store> &dataset;  ///< dataset to write moments to
  WriteGbx writer;          ///< pointer to attribute ragged arrays in dataset

  // use functor from writer to collect data from gridboxes in parallel
  void collect_data_from_gridboxes(const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    auto functor = writer.get_functor(d_gbxs);
    Kokkos::parallel_for("collect_superdrops_data", Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
                         functor);
    // TODO(CB) ragged count
  }

  // collect superdroplet data from gridboxes and then write it to arrays in the dataset
  void at_start_step(const viewd_constgbx d_gbxs) const {
    collect_data_from_gridboxes(d_gbxs);
    writer.write_to_array(dataset);
  }

 public:
  DoSuperdropsObs(Dataset<Store> &dataset, WriteGbx writer) : dataset(dataset), writer(writer) {}

  ~DoSuperdropsObs() { writer.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes superdrops observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs) const {
    at_start_step(d_gbxs);
  }
};

/* constructs observer which writes mass moments of droplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto SuperdropsObserver(const unsigned int interval, Dataset<Store> &dataset,
                                        const int maxchunk, const size_t ngbxs) {
  return ConstTstepObserver(interval, DoSuperdropsObs(dataset, writer));
}

#endif  // LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_
