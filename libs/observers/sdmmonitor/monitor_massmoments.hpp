/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_massmoments.hpp
 * Project: sdmmonitor
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 27th May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct to create observer which outputs the average mass moments monitored
 * from SDM microphysical process in each gridbox a constant interval at the
 * start of each timestep.
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>

#include "../../kokkosaliases.hpp"
#include "../consttstep_observer.hpp"
#include "../observers.hpp"
#include "./do_sdmmonitor_obs.hpp"

/* struct satisfies SDMMonitor concept for use in do_sdmmonitor_obs to make observer */
struct MonitorMassMoments {
  using datatype = float;
  size_t monitor_count;  // number of calls to monitor microphysics since last reset
  Buffer<datatype>::mirrorviewd_buffer d_data;  // view on device copied to host by DoSDMMonitorObs

  /**
   * @brief Parallel loop to fill d_data with zero value.
   */
  void reset_monitor() const;

  /**
   * @brief Placeholder function to obey SDMMonitor concept does nothing.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param totmass_condensed Mass condensed in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember& team_member, const double totmass_condensed) const {}

  /**
   * @brief Monitor 0th, 1st and 2nd moments of the droplet mass distribution
   *
   * Averages current of mass moments with value stored since d_data was last reset.
   *
   * _Note:_ possible conversion of mass moments at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output depending on datatype alias.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember& team_member, const viewd_constsupers supers) const;

  /**
   * @brief Constructor for MonitorMassMoments
   *
   * @param ngbxs Number of gridboxes in domain.
   */
  explicit MonitorMassMoments(const size_t ngbxs) : d_data("massmom_todo", ngbxs) {
    reset_monitor();
  }
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
inline Observer auto MonitorMassMoments(const unsigned int interval, Dataset<Store>& dataset,
                                        const size_t maxchunk, const size_t ngbxs) {
  using Mo = MonitorMassMoments;
  const auto name = std::string_view("massmom_todo");
  const auto units = std::string_view("todo");
  constexpr auto scale_factor = 1.0;  // TODO(CB): appropriate metadata
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr_ptr = std::make_shared<XarrayZarrArray<Store, Mo::datatype>>(
      dataset.template create_array<Mo::datatype>(name, units, scale_factor, chunkshape, dimnames));

  const auto do_obs = DoSDMMonitorObs<Store, Mo, Mo::datatype>(dataset, xzarr_ptr, Mo(ngbxs));
  return ConstTstepObserver(interval, do_obs);
}

#endif  //  LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_HPP_
