/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: state_observer.hpp
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

#ifndef LIBS_OBSERVERS2_STATE_OBSERVER_HPP_
#define LIBS_OBSERVERS2_STATE_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <cassert>
#include <concepts>
#include <iostream>
#include <memory>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "./write_gridboxes.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"

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

  void write_to_array(Dataset<Store> &dataset) const {
    Kokkos::deep_copy(h_data, d_data);
    dataset.write_to_array(xzarr_ptr, h_data);
  }

  void write_arrayshape(Dataset<Store> &dataset) const { dataset.write_arrayshape(xzarr_ptr); }
};

/* constructs observer which writes variables from the state of each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto StateObserver(const unsigned int interval, Dataset<Store> &dataset,
                                   const int maxchunk, const size_t ngbxs) {
  const auto write_state = DataFromGridboxesToArray(dataset, maxchunk, ngbxs);
  return ConstTstepObserver(interval, WriteGridboxes(dataset, write_state));
}

#endif  // LIBS_OBSERVERS2_STATE_OBSERVER_HPP_
