/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: time_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output time at the start of each observation timestep as a coordinate of an Xarray
 * dataset
 */

#ifndef LIBS_OBSERVERS_TIME_OBSERVER_HPP_
#define LIBS_OBSERVERS_TIME_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/observers.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "zarr/xarray_zarr_array.hpp"

/**
 * @class DoTimeObs
 * @brief Template class for functionality to observe time at the start of each timestep and write
 * it to a Zarr array as a coordinate of an Xarray dataset.
 * @tparam Dataset Type of dataset.
 * @tparam Store Type of store which dataset writes to.
 */
template <typename Dataset, typename Store>
class DoTimeObs {
 private:
  Dataset &dataset; /**< Dataset to write time data to. */
  std::shared_ptr<XarrayZarrArray<Store, float>>
      xzarr_ptr; /**< Pointer to time array in dataset. */
  std::function<double(unsigned int)>
      step2dimlesstime; /**< Function to convert timesteps to real time [assumed seconds]. */

  /**
   * @brief Increment size of time dimension in dataset and write out current time of model (assumed
   * seconds) to the array in the dataset.
   *
   * _Note:_ conversion of time from double precision (8 bytes double) to single precision (4 bytes
   * float) in output.
   *
   * @param t_mdl Current model time.
   */
  void at_start_step(const unsigned int t_mdl) const {
    const auto ntimes = size_t{dataset.get_dimension("time") + 1};
    const auto timedim = std::pair<std::string, size_t>({"time", ntimes});
    dataset.set_dimension(timedim);

    const auto time = static_cast<float>(step2dimlesstime(t_mdl));
    dataset.write_to_array(xzarr_ptr, time);
  }

 public:
  /**
   * @brief Constructor for DoTimeObs.
   * @param dataset Dataset to write time data to.
   * @param store Store which dataset writes to.
   * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
   * @param step2dimlesstime Function to convert model timesteps to a real time [assumed seconds].
   */
  DoTimeObs(Dataset &dataset, Store &store, const size_t maxchunk,
            const std::function<double(unsigned int)> step2dimlesstime)
      : dataset(dataset),
        xzarr_ptr(std::make_shared<XarrayZarrArray<Store, float>>(
            dataset.template create_coordinate_array<float>("time", "s", dlc::TIME0, maxchunk, 0))),
        step2dimlesstime(step2dimlesstime) {}

  /**
   * @brief Destructor for DoTimeObs.
   */
  ~DoTimeObs() { dataset.write_arrayshape(xzarr_ptr); }

  /**
   * @brief Placeholder for before timestepping functionality and to make class satisfy observer
   * concept.
   */
  void before_timestepping(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
    std::cout << "observer includes time observer\n";
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
   * @param d_supers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const subviewd_constsupers d_supers) const {
    at_start_step(t_mdl);
  }

  /**
   * @brief Get null monitor for SDM processes from observer.
   *
   * @return monitor 'mo' of the observer that does nothing
   */
  SDMMonitor auto get_sdmmonitor() const { return NullSDMMonitor{}; }
};

/**
 * @brief Constructs an observer which writes (real) time at start of each observation timestep to a
 * 1-D array with a constant observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @tparam Dataset Type of dataset
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param store Store which dataset writes to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param step2dimlesstime Function to convert model timesteps to real time (assumed seconds).
 * @return Constructed type satisfying observer concept.
 */
template <typename Dataset, typename Store>
inline Observer auto TimeObserver(const unsigned int interval, Dataset &dataset, Store &store,
                                  const size_t maxchunk,
                                  const std::function<double(unsigned int)> step2dimlesstime) {
  return ConstTstepObserver(interval, DoTimeObs(dataset, store, maxchunk, step2dimlesstime));
}

#endif  // LIBS_OBSERVERS_TIME_OBSERVER_HPP_
