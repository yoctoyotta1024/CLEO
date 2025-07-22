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
#include "zarr/xarray_zarr_array.hpp"

namespace KCS = KokkosCleoSettings;

/**
 * @class DoSDMMonitorObs
 * @brief Template class for functionality to observe data from a monitor of a SDM process at the
 * start of each timestep and write it to a Zarr array in an Xarray dataset.
 * @tparam Store Type of store for dataset.
 */
template <typename Dataset, typename Store, SDMMonitor SDMMo, typename T>
class DoSDMMonitorObs {
 private:
  using viewh_buffer = Buffer<T>::viewh_buffer;
  Dataset &dataset;                                     /**< Dataset to write time data to. */
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
  DoSDMMonitorObs(Dataset &dataset, Store &store,
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

#endif  // LIBS_OBSERVERS_SDMMONITOR_DO_SDMMONITOR_OBS_HPP_
