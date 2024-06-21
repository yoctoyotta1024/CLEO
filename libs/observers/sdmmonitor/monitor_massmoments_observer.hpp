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
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output variables from a mass moments monitor of SDM processes at
 * a constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_OBSERVER_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../../kokkosaliases.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/create_massmoments_arrays.hpp"
#include "observers/observers.hpp"
#include "observers/sdmmonitor/monitor_massmoments.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "zarr/buffer.hpp"
#include "zarr/dataset.hpp"
#include "zarr/xarray_zarr_array.hpp"

template <typename Store>
struct MonitorMassMomentXarrays {
  XarrayZarrArray<Store, uint64_t> mom0_microphys; /**< 0th mass moment from microphysics Xarray */
  XarrayZarrArray<Store, float> mom1_microphys;    /**< 1st mass moment from microphysics Xarray */
  XarrayZarrArray<Store, float> mom2_microphys;    /**< 2nd mass moment from microphysics Xarray */
  XarrayZarrArray<Store, uint64_t> mom0_motion;    /**< 0th mass moment from motion Xarray */
  XarrayZarrArray<Store, float> mom1_motion;       /**< 1st mass moment from motion Xarray */
  XarrayZarrArray<Store, float> mom2_motion;       /**< 2nd mass moment from motion Xarray */

  MonitorMassMomentXarrays(const Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs)
      : mom0_microphys(create_massmom0_xarray(dataset, "massmom0_microphys", maxchunk, ngbxs)),
        mom1_microphys(create_massmom1_xarray(dataset, "massmom1_microphys", maxchunk, ngbxs)),
        mom2_microphys(create_massmom2_xarray(dataset, "massmom2_microphys", maxchunk, ngbxs)),
        mom0_motion(create_massmom0_xarray(dataset, "massmom0_motion", maxchunk, ngbxs)),
        mom1_motion(create_massmom1_xarray(dataset, "massmom1_motion", maxchunk, ngbxs)),
        mom2_motion(create_massmom2_xarray(dataset, "massmom2_motion", maxchunk, ngbxs)) {}
};

/**
 * @class DoMonitorMassMomentsObs
 * @brief Class for functionality to observe data from a mass moments monitor of a SDM
 * process at the start of each timestep and write it to a Zarr array in an Xarray dataset.
 * @tparam Store Type of store for dataset.
 */
template <typename Store, typename MonitorViewsType>
class DoMonitorMassMomentsObs {
 private:
  Dataset<Store> &dataset; /**< Dataset to write time data to. */
  std::shared_ptr<MonitorMassMomentXarrays<Store>> xzarrs_ptr; /**< Pointer to arrays in dataset. */
  MonitorMassMoments<MonitorViewsType> monitor;

  /**
   * @brief Copy data from d_data view on device into host view,
   * then write to the array in the dataset.
   */
  template <typename T>
  void write_to_array(const Buffer<T>::mirrorviewd_buffer d_data,
                      XarrayZarrArray<Store, T> &xzarr) const {
    using viewh_buffer = Buffer<T>::viewh_buffer;
    const auto h_data = viewh_buffer("h_data", d_data.extent(0));
    Kokkos::deep_copy(h_data, d_data);
    dataset.write_to_array(xzarr, h_data);
  }

  /**
   * @brief Write each mass moment from the monitor's views to the appropriate arrays in the dataset
   * then reset the monitor.
   */
  void at_start_step() const {
    write_to_array(monitor.microphysics_moms.d_mom0, xzarrs_ptr->mom0_microphys);
    write_to_array(monitor.microphysics_moms.d_mom1, xzarrs_ptr->mom1_microphys);
    write_to_array(monitor.microphysics_moms.d_mom2, xzarrs_ptr->mom2_microphys);
    write_to_array(monitor.motion_moms.d_mom0, xzarrs_ptr->mom0_motion);
    write_to_array(monitor.motion_moms.d_mom1, xzarrs_ptr->mom1_motion);
    write_to_array(monitor.motion_moms.d_mom2, xzarrs_ptr->mom2_motion);

    monitor.reset_monitor();
  }

 public:
  /**
   * @brief Constructor for DoMonitorMassMomentsObs.
   * @param dataset Dataset to write monitored data to.
   * @param maxchunk The maximum chunk size (number of elements) for Xarrays.
   * @param ngbxs The number of gridboxes.
   */
  DoMonitorMassMomentsObs(Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs)
      : dataset(dataset),
        xzarrs_ptr(std::make_shared<MonitorMassMomentXarrays<Store>>(dataset, maxchunk, ngbxs)),
        monitor(ngbxs) {}

  /**
   * @brief Destructor for DoMonitorMassMomentsObs.
   */
  ~DoMonitorMassMomentsObs() {
    dataset.write_arrayshape(xzarrs_ptr->mom0_microphys);
    dataset.write_arrayshape(xzarrs_ptr->mom1_microphys);
    dataset.write_arrayshape(xzarrs_ptr->mom2_microphys);
    dataset.write_arrayshape(xzarrs_ptr->mom0_motion);
    dataset.write_arrayshape(xzarrs_ptr->mom1_motion);
    dataset.write_arrayshape(xzarrs_ptr->mom2_motion);
  }

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
 * @param ngbxs The number of gridboxes.
 * @return Constructed type satisfying observer concept.
 */
template <typename Store>
inline Observer auto MonitorMassMomentsObserver(const unsigned int interval,
                                                Dataset<Store> &dataset, const size_t maxchunk,
                                                const size_t ngbxs) {
  const auto do_obs =
      DoMonitorMassMomentsObs<Store, MonitorMassMomentViews>(dataset, maxchunk, ngbxs);
  return ConstTstepObserver(interval, do_obs);
}

#endif  // LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_OBSERVER_HPP_
