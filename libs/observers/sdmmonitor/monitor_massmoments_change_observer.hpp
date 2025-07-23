/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_massmoments_change_observer.hpp
 * Project: sdmmonitor
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output variables from a mass moments monitor of SDM processes at
 * a constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_CHANGE_OBSERVER_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_CHANGE_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <memory>

#include "../../kokkosaliases.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/create_massmoments_arrays.hpp"
#include "observers/observers.hpp"
#include "observers/sdmmonitor/monitor_massmoments_change.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "zarr/buffer.hpp"
#include "zarr/xarray_zarr_array.hpp"

namespace KCS = KokkosCleoSettings;

template <typename Dataset, typename Store>
struct MonitorMassMomentsChangeXarrays {
  XarrayZarrArray<Store, uint64_t>
      delta_mom0_microphys; /**< change in 0th mass moment from microphysics Xarray */
  XarrayZarrArray<Store, float>
      delta_mom1_microphys; /**< change in 1st mass moment from microphysics Xarray */
  XarrayZarrArray<Store, float>
      delta_mom2_microphys; /**< change in 2nd mass moment from microphysics Xarray */
  XarrayZarrArray<Store, uint64_t>
      delta_mom0_motion; /**< change in 0th mass moment from motion Xarray */
  XarrayZarrArray<Store, float>
      delta_mom1_motion; /**< change in 1st mass moment from motion Xarray */
  XarrayZarrArray<Store, float>
      delta_mom2_motion; /**< change in 2nd mass moment from motion Xarray */

  MonitorMassMomentsChangeXarrays(const Dataset &dataset, Store &store, const size_t maxchunk,
                                  const size_t ngbxs)
      : delta_mom0_microphys(
            create_massmom0_xarray(dataset, store, "delta_massmom0_microphys", maxchunk, ngbxs)),
        delta_mom1_microphys(
            create_massmom1_xarray(dataset, store, "delta_massmom1_microphys", maxchunk, ngbxs)),
        delta_mom2_microphys(
            create_massmom2_xarray(dataset, store, "delta_massmom2_microphys", maxchunk, ngbxs)),
        delta_mom0_motion(
            create_massmom0_xarray(dataset, store, "delta_massmom0_motion", maxchunk, ngbxs)),
        delta_mom1_motion(
            create_massmom1_xarray(dataset, store, "delta_massmom1_motion", maxchunk, ngbxs)),
        delta_mom2_motion(
            create_massmom2_xarray(dataset, store, "delta_massmom2_motion", maxchunk, ngbxs)) {}
};

template <typename Dataset, typename Store>
struct MonitorRainMassMomentsChangeXarrays {
  XarrayZarrArray<Store, uint64_t>
      delta_mom0_microphys; /**< change in 0th mass moment from microphysics Xarray */
  XarrayZarrArray<Store, float>
      delta_mom1_microphys; /**< change in 1st mass moment from microphysics Xarray */
  XarrayZarrArray<Store, float>
      delta_mom2_microphys; /**< change in 2nd mass moment from microphysics Xarray */
  XarrayZarrArray<Store, uint64_t>
      delta_mom0_motion; /**< change in 0th mass moment from motion Xarray */
  XarrayZarrArray<Store, float>
      delta_mom1_motion; /**< change in 1st mass moment from motion Xarray */
  XarrayZarrArray<Store, float>
      delta_mom2_motion; /**< change in 2nd mass moment from motion Xarray */

  MonitorRainMassMomentsChangeXarrays(const Dataset &dataset, Store &store, const size_t maxchunk,
                                      const size_t ngbxs)
      : delta_mom0_microphys(create_massmom0_xarray(
            dataset, store, "delta_massmom0_raindrops_microphys", maxchunk, ngbxs)),
        delta_mom1_microphys(create_massmom1_xarray(
            dataset, store, "delta_massmom1_raindrops_microphys", maxchunk, ngbxs)),
        delta_mom2_microphys(create_massmom2_xarray(
            dataset, store, "delta_massmom2_raindrops_microphys", maxchunk, ngbxs)),
        delta_mom0_motion(create_massmom0_xarray(dataset, store, "delta_massmom0_raindrops_motion",
                                                 maxchunk, ngbxs)),
        delta_mom1_motion(create_massmom1_xarray(dataset, store, "delta_massmom1_raindrops_motion",
                                                 maxchunk, ngbxs)),
        delta_mom2_motion(create_massmom2_xarray(dataset, store, "delta_massmom2_raindrops_motion",
                                                 maxchunk, ngbxs)) {}
};

/**
 * @class DoMonitorMassMomentsChangeObs
 * @brief Class for functionality to observe data from a mass moments monitor of a SDM
 * process at the start of each timestep and write it to a Zarr array in an Xarray dataset.
 * @tparam Store Type of store for dataset.
 * @tparam MonitorXarraysType Type for xarrays for mass moments in dataset.
 * @tparam MonitorViewsType Type for views which calculate mass moments for xarrays.
 */
template <typename Dataset, typename Store, typename MonitorXarraysType, typename MonitorViewsType>
class DoMonitorMassMomentsChangeObs {
 private:
  Dataset &dataset;                               /**< Dataset to write time data to. */
  std::shared_ptr<MonitorXarraysType> xzarrs_ptr; /**< Pointer to arrays in dataset. */
  MonitorMassMomentsChange<MonitorViewsType> monitor;

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
   * @brief Write each change in mass moment from the monitor's views to the
   * appropriate arrays in the dataset then reset the monitor.
   */
  void at_start_step() const {
    write_to_array(monitor.microphysics_moms.d_delta_mom0, xzarrs_ptr->delta_mom0_microphys);
    write_to_array(monitor.microphysics_moms.d_delta_mom1, xzarrs_ptr->delta_mom1_microphys);
    write_to_array(monitor.microphysics_moms.d_delta_mom2, xzarrs_ptr->delta_mom2_microphys);
    write_to_array(monitor.motion_moms.d_delta_mom0, xzarrs_ptr->delta_mom0_motion);
    write_to_array(monitor.motion_moms.d_delta_mom1, xzarrs_ptr->delta_mom1_motion);
    write_to_array(monitor.motion_moms.d_delta_mom2, xzarrs_ptr->delta_mom2_motion);

    monitor.reset_monitor();
  }

 public:
  /**
   * @brief Constructor for DoMonitorMassMomentsChangeObs.
   * @param dataset Dataset to write monitored data to.
   * &param store Store dataset writes into.
   * @param maxchunk The maximum chunk size (number of elements) for Xarrays.
   * @param ngbxs The number of gridboxes.
   */
  DoMonitorMassMomentsChangeObs(Dataset &dataset, Store &store, const size_t maxchunk,
                                const size_t ngbxs)
      : dataset(dataset),
        xzarrs_ptr(std::make_shared<MonitorXarraysType>(dataset, store, maxchunk, ngbxs)),
        monitor(ngbxs) {}

  /**
   * @brief Destructor for DoMonitorMassMomentsChangeObs.
   */
  ~DoMonitorMassMomentsChangeObs() {
    dataset.write_arrayshape(xzarrs_ptr->delta_mom0_microphys);
    dataset.write_arrayshape(xzarrs_ptr->delta_mom1_microphys);
    dataset.write_arrayshape(xzarrs_ptr->delta_mom2_microphys);
    dataset.write_arrayshape(xzarrs_ptr->delta_mom0_motion);
    dataset.write_arrayshape(xzarrs_ptr->delta_mom1_motion);
    dataset.write_arrayshape(xzarrs_ptr->delta_mom2_motion);
  }

  /**
   * @brief Placeholder for before timestepping functionality and to make class satisfy observer
   * concept.
   */
  void before_timestepping(const viewd_constgbx d_gbxs,
                           const subviewd_constsupers domainsupers) const {
    std::cout << "observer includes SDM monitor observer\n";

    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for(
        "monitor_before_timestepping", TeamPolicy(ngbxs, KCS::team_size),
        KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
          const auto ii = team_member.league_rank();
          const auto supers = d_gbxs(ii).supersingbx.readonly(domainsupers);
          monitor.before_timestepping(team_member, supers);
        });
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
   * @param d_supers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const subviewd_constsupers d_supers) const {
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
 * @param store Store which dataset writes to
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Constructed type satisfying observer concept.
 */
template <typename Dataset, typename Store>
inline Observer auto MonitorMassMomentsChangeObserver(const unsigned int interval, Dataset &dataset,
                                                      Store &store, const size_t maxchunk,
                                                      const size_t ngbxs) {
  const auto do_obs =
      DoMonitorMassMomentsChangeObs<Dataset, Store, MonitorMassMomentsChangeXarrays<Dataset, Store>,
                                    MonitorMassMomentsChangeViews>(dataset, store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, do_obs);
}

/**
 * @brief Constructs an observer which writes data monitoring the mass moments of the raindrop's
 * distributions during microphysics and super-droplet motion to arrays with a constant observation
 * timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param store Store which dataset writes to
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Constructed type satisfying observer concept.
 */
template <typename Dataset, typename Store>
inline Observer auto MonitorRainMassMomentsObserver(const unsigned int interval, Dataset &dataset,
                                                    Store &store, const size_t maxchunk,
                                                    const size_t ngbxs) {
  const auto do_obs =
      DoMonitorMassMomentsChangeObs<Dataset, Store,
                                    MonitorRainMassMomentsChangeXarrays<Dataset, Store>,
                                    MonitorRainMassMomentsChangeViews>(dataset, store, maxchunk,
                                                                       ngbxs);
  return ConstTstepObserver(interval, do_obs);
}

#endif  // LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_CHANGE_OBSERVER_HPP_
