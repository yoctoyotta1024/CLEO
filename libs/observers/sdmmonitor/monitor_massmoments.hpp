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
 * Last Modified: Tuesday 28th May 2024
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

/* struct satisfies SDMMonitor concept for use in do_sdmmonitor_obs to make observer */
struct MonitorMassMoments {
  using datatype = float;
  using viewd_count = Kokkos::View<size_t[1]>;
  viewd_count microphysics_count;  // number of calls to monitor microphysics since last reset
  viewd_count motion_count;        // number of calls to monitor motion since last reset
  Buffer<datatype>::mirrorviewd_buffer d_data;  // view on device copied to host by DoSDMMonitorObs

  /**
   * @brief Parallel loop to fill d_data with zero value and set monitor_XXX_counts to zero.
   */
  void reset_monitor() const;

  /**
   * @brief Write the 0th, 1st and 2nd moments of the droplet mass distribution to data views.
   *
   * Calculates the current mass moments and then averages them with the current values for the
   * mass moments stored since the data views were last reset.
   *
   * _Note:_ possible conversion of mass moments at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output depending on datatype alias.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  size_t average_massmoments(const TeamMember& team_member, const viewd_constsupers supers,
                             size_t count) const;

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
   * calls average_massmoments to monitor the moments of the droplet mass
   * distribution during SDM microphysics
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember& team_member, const viewd_constsupers supers) const {
    microphysics_count(0) = average_massmoments(team_member, supers, microphysics_count(0));
  }

  /**
   * @brief Monitor 0th, 1st and 2nd moments of the droplet mass distribution
   *
   * calls average_massmoments to monitor the moments of the droplet mass
   * distribution during SDM motion.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one motion timestep
   */
  KOKKOS_FUNCTION
  void monitor_motion(const TeamMember& team_member, const viewd_constsupers supers) const {
    motion_count(0) = average_massmoments(team_member, supers, motion_count(0));
  }

  /**
   * @brief Constructor for MonitorMassMoments
   *
   * @param ngbxs Number of gridboxes in domain.
   */
  explicit MonitorMassMoments(const size_t ngbxs) : d_data("massmom_todo", ngbxs) {
    reset_monitor();
  }
};

#endif  //  LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_HPP_
