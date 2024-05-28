/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_massmoments_observer.hpp
 * Project: sdmmonitor
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 28th May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output variables from a mass moments monitor of SDM processes at a constant interval
 * at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_OBSERVER_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../../kokkosaliases.hpp"
#include "../consttstep_observer.hpp"
#include "../observers.hpp"
#include "./monitor_massmoments.hpp"
#include "zarr/buffer.hpp"
#include "zarr/dataset.hpp"
#include "zarr/xarray_zarr_array.hpp"

/**
 * @class DoMonitorMassMomentsObs
 * @brief Class for functionality to observe data from a mass moments monitor of a SDM
 * process at the start of each timestep and write it to a Zarr array in an Xarray dataset.
 * @tparam Store Type of store for dataset.
 */
template <typename Store, typename T>
class DoMonitorMassMomentsObs {
 private:
  using viewh_buffer = Buffer<T>::viewh_buffer;
  Dataset<Store> &dataset;                              /**< Dataset to write time data to. */
  std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr; /**< Pointer to array in dataset. */
  MonitorMassMoments monitor;

  /**
   * @brief Copy data from monitor to the array in the dataset then reset monitor.
   */
  void at_start_step() const {
    const auto h_data = viewh_buffer("h_data", monitor.d_data.extent(0));
    Kokkos::deep_copy(h_data, monitor.d_data);
    dataset.write_to_array(xzarr_ptr, h_data);

    monitor.reset_monitor();
  }

 public:
  /**
   * @brief Constructor for DoMonitorMassMomentsObs.
   * @param dataset Dataset to write monitored data to.
   * @param xzarr_ptr Pointer to zarr array in xarray dataset.
   */
  DoMonitorMassMomentsObs(Dataset<Store> &dataset,
                          const std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr,
                          const MonitorMassMoments monitor)
      : dataset(dataset), xzarr_ptr(xzarr_ptr), monitor(monitor) {}

  /**
   * @brief Destructor for DoMonitorMassMomentsObs.
   */
  ~DoMonitorMassMomentsObs() { dataset.write_arrayshape(xzarr_ptr); }

  /**
   * @brief Placeholder for before timestepping functionality and to make class satisfy observer
   * concept.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes SDM monitor observer\n";
  }

  /**
   * @brief Placeholder for after timestepping functionality and to make class satisfy observer
   * concept.
   */
  void after_timestepping() const {}

  /**
   * @brief Adapter to call at start step function which writes data from the monitor to the array
   * in the dataset.
   *
   * @param t_mdl Current model timestep.
   * @param d_gbxs View of gridboxes on device.
   * @param totsupers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step();
  }

  /**
   * @brief Get monitor for SDM processes from observer.
   *
   * @return monitor 'mo' of the observer
   */
  SDMMonitor auto get_sdmmonitor() const { return monitor; }
};

/**
 * @brief Constructs an observer which writes data monitoring the mass moments during microphysics
 * and super-droplet motion to arrays with a constant observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @return Constructed type satisfying observer concept.
 */
template <typename Store>
inline Observer auto MonitorMassMomentsObserver(const unsigned int interval,
                                                Dataset<Store> &dataset, const size_t maxchunk,
                                                const size_t ngbxs) {
  using Mo = MonitorMassMoments;
  const auto name = std::string_view("massmom_todo");
  const auto units = std::string_view("todo");
  constexpr auto scale_factor = 1.0;  // TODO(CB): appropriate metadata
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr_ptr = std::make_shared<XarrayZarrArray<Store, Mo::datatype>>(
      dataset.template create_array<Mo::datatype>(name, units, scale_factor, chunkshape, dimnames));

  const auto do_obs =
      DoMonitorMassMomentsObs<Store, Mo, Mo::datatype>(dataset, xzarr_ptr, Mo(ngbxs));
  return ConstTstepObserver(interval, do_obs);
}

#endif  // LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_OBSERVER_HPP_
