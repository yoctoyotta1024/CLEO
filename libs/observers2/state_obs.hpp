/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: state_obs.hpp
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
 * Observer to output variables related to Gridboxes' state at the start of
 * each timestep to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_STATE_OBS_HPP_
#define LIBS_OBSERVERS2_STATE_OBS_HPP_

#include <Kokkos_Core.hpp>
#include <cassert>
#include <concepts>
#include <iostream>
#include <memory>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"
#include "zarr2/dataset.hpp"

// Functor to perform copy in parallel of 1 value (pressure) from each gridbox
template <typename Store>
class DataFromGridboxesToArray {
 private:
  using viewh_data = Kokkos::View<double *, HostSpace>;
  using mirrorviewd_data = Kokkos::View<double *, HostSpace::array_layout,
                                        ExecSpace>;  // type for mirror of host view on device
  std::shared_ptr<XarrayZarrArray<Store, double>> xzarr_ptr;
  viewh_data h_data;
  mirrorviewd_data d_data;

 public:
  struct Functor {
    viewd_constgbx d_gbxs;
    mirrorviewd_data d_data;

    Functor(const viewd_constgbx d_gbxs, mirrorviewd_data d_data)
        : d_gbxs(d_gbxs), d_data(d_data) {}

    // Functor operator to perform copy of each element in parallel
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const { d_data(ii) = d_gbxs(ii).state.press; }
  };

  // Constructor to initialize Kokkos view
  DataFromGridboxesToArray(Dataset<Store> &dataset, const int maxchunk, const size_t ngbxs)
      : xzarr_ptr(
            std::make_shared<XarrayZarrArray<Store, double>>(dataset.template create_array<double>(
                "press", "hPa", "<f8", dlc::P0 / 100, good2Dchunkshape(maxchunk, ngbxs),
                {"time", "gbxindex"}))),
        h_data("h_data", ngbxs),
        d_data(Kokkos::create_mirror_view(ExecSpace(), h_data)) {}

  Functor get_functor(const viewd_constgbx d_gbxs) const {
    assert((d_gbxs.extent(0) == d_data.extent(0)) &&
           "d_data view must be size of the number of gridboxes");
    return Functor(d_gbxs, d_data);
  }

  void write_data(Dataset<Store> &dataset) const {
    Kokkos::deep_copy(h_data, d_data);
    dataset.write_to_array(xzarr_ptr, h_data);
  }

  void write_arrayshape(Dataset<Store> &dataset) const { dataset.write_arrayshape(xzarr_ptr); }
};

/* observe variables in the state of each
gridbox and write them to repspective arrays
in a store as determined by the Dataset */
template <typename Store>
class DoStateObs {
 private:
  Dataset<Store> &dataset;
  DataFromGridboxesToArray<Store> data2array;

  void fetch_data_from_gridboxes(const viewd_constgbx d_gbxs) const {
    auto functor = data2array.get_functor(d_gbxs);

    const auto ngbxs = size_t{d_gbxs.extent(0)};
    Kokkos::parallel_for("stateobs", Kokkos::RangePolicy<ExecSpace>(0, ngbxs), functor);
  }

 public:
  DoStateObs(Dataset<Store> &dataset, const int maxchunk, const size_t ngbxs)
      : dataset(dataset), data2array(dataset, maxchunk, ngbxs) {}

  ~DoStateObs() { data2array.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes State observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    fetch_data_from_gridboxes(d_gbxs);
    data2array.write_data(dataset);

    // dataset.set_dimension({"time", time+1}); // TODO(CB) do this with coord observer
  }
};

/* constructs observer of variables in the state
of each gridbox with a constant timestep 'interval'
using an instance of the DoStateObs class */
template <typename Store>
inline Observer auto StateObserver(const unsigned int interval, Dataset<Store> &dataset,
                                   const int maxchunk, const size_t ngbxs) {
  const auto obs = DoStateObs<Store>(dataset, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

#endif  // LIBS_OBSERVERS2_STATE_OBS_HPP_
