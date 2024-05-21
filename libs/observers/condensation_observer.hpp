/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: condensation_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 21st May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output condensation rate monitored from SDM microphysical process in each gridbox a
 * constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS_CONDENSATION_OBSERVER_HPP_
#define LIBS_OBSERVERS_CONDENSATION_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "../kokkosaliases.hpp"
#include "./consttstep_observer.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "sdmmonitor/monitor_condensation.hpp"
#include "sdmmonitor/sdmmonitor.hpp"
#include "zarr/dataset.hpp"
#include "zarr/xarray_zarr_array.hpp"

/**
 * @class DoCondensationObs
 * @brief Template class for functionality to observe monitor of condensation microphysics at the
 * start of each timestep and write it to a Zarr array in an Xarray dataset.
 * @tparam Store Type of store for dataset.
 */
template <typename Store>
class DoCondensationObs {
 private:
  Dataset<Store> &dataset; /**< Dataset to write time data to. */
  std::shared_ptr<XarrayZarrArray<Store, float>>
      xzarr_ptr; /**< Pointer to condrate array in dataset. */
  MonitorCondensation monitor;

  /**
   * @brief Copy data from monitor to the array in the dataset.
   */
  void at_start_step() const {
    const size_t sz = monitor.condrate.extent(0);
    auto h_data = Buffer<float>::viewh_buffer("h_data", sz);
    Kokkos::deep_copy(h_data, monitor.condrate);
    dataset.write_to_array(xzarr_ptr, h_data);
  }

 public:
  /**
   * @brief Constructor for DoCondensationObs.
   * @param dataset Dataset to write condensation data to.
   * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
   */
  DoCondensationObs(Dataset<Store> &dataset, const size_t maxchunk)
      : dataset(dataset),
        xzarr_ptr(
            std::make_shared<XarrayZarrArray<Store, float>>(dataset.template create_array<float>(
                "condrate", "TODO(CB)", "<f4", 0.5, maxchunk, "time"))) {}

  /**
   * @brief Destructor for DoCondensationObs.
   */
  ~DoCondensationObs() { dataset.write_arrayshape(xzarr_ptr); }

  /**
   * @brief Placeholder for before timestepping functionality and to make class satisfy observer
   * concept.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes condensation observer\n";
  }

  /**
   * @brief Placeholder for after timestepping functionality and to make class satisfy observer
   * concept.
   */
  void after_timestepping() const {}

  /**
   * @brief Adapter to call at start step function which writes the current time of the model
   * (assumed seconds) to the array in the dataset.
   *
   * @param t_mdl Current model timestep.
   * @param d_gbxs View of gridboxes on device.
   * @param totsupers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step();
  }

  SDMMonitor auto get_sdmmonitor() const { return monitor; }
};

/**
 * @brief Constructs an observer which writes data monitoring condensation microphysics to an
 * array with a constant observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @return Constructed type satisfying observer concept.
 */
template <typename Store>
inline Observer auto CondensationObserver(const unsigned int interval, Dataset<Store> &dataset,
                                          const size_t maxchunk) {
  return ConstTstepObserver(interval, DoCondensationObs(dataset, maxchunk));
}

#endif  // LIBS_OBSERVERS_CONDENSATION_OBSERVER_HPP_
