/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: do_sdmmonitor_obs.hpp
 * Project: sdmmonitor
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 15th June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output variables from a monitor of SDM processes at a constant interval at the start
 * of each timestep.
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_DO_SDMMONITOR_OBS_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_DO_SDMMONITOR_OBS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../../kokkosaliases.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "zarr/buffer.hpp"
#include "zarr/dataset.hpp"
#include "zarr/xarray_zarr_array.hpp"

/**
 * @class DoSDMMonitorObs
 * @brief Template class for functionality to observe data from a monitor of a SDM process at the
 * start of each timestep and write it to a Zarr array in an Xarray dataset.
 * @tparam Store Type of store for dataset.
 */
template <typename Store, SDMMonitor SDMMo, typename T>
class DoSDMMonitorObs {
 private:
  using viewh_buffer = Buffer<T>::viewh_buffer;
  Dataset<Store> &dataset;                              /**< Dataset to write time data to. */
  std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr; /**< Pointer to array in dataset. */
  SDMMo monitor;

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
   * @brief Constructor for DoSDMMonitorObs.
   * @param dataset Dataset to write monitored data to.
   * @param xzarr_ptr Pointer to zarr array in xarray dataset.
   * @param monitor SDMMonitor to use.
   */
  DoSDMMonitorObs(Dataset<Store> &dataset,
                  const std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr, const SDMMo monitor)
      : dataset(dataset), xzarr_ptr(xzarr_ptr), monitor(monitor) {}

  /**
   * @brief Destructor for DoSDMMonitorObs.
   */
  ~DoSDMMonitorObs() { dataset.write_arrayshape(xzarr_ptr); }

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

#endif  // LIBS_OBSERVERS_SDMMONITOR_DO_SDMMONITOR_OBS_HPP_
