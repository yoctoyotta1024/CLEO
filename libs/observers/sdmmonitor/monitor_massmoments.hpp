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
 * Last Modified: Thursday 6th June 2024
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
#include "zarr/buffer.hpp"

struct MonitorMassMomentViews {
  Buffer<uint64_t>::mirrorviewd_buffer d_massmom0;  // view on device for monitoring 0th mass moment
  Buffer<float>::mirrorviewd_buffer d_massmom1;     // view on device for monitoring 1st mass moment
  Buffer<float>::mirrorviewd_buffer d_massmom2;     // view on device for monitoring 2nd mass moment

  /**
   * @brief Parallel loop to fill device views with zero value
   */
  void reset_views() const;

  /**
   * @brief Write the 0th, 1st and 2nd moments of the droplet mass distribution to data views.
   *
   * Calculates the current mass moments and then overwrites the current values for the
   * mass moments stored since the data views were last reset.
   *
   * _Note:_ possible conversion of mass moments at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox
   */
  KOKKOS_FUNCTION
  void calculate_massmoments(const TeamMember& team_member, const viewd_constsupers supers) const;

  explicit MonitorMassMomentViews(const size_t ngbxs)
      : d_massmom0("d_monitor_massmom0", ngbxs),
        d_massmom1("d_monitor_massmom1", ngbxs),
        d_massmom2("d_monitor_massmom2", ngbxs) {
    reset_views();
  }
};

/* struct satisfies SDMMonitor concept in order to make observer for monitoring mass moments
 * according to the templated MonitorViewsType e.g. 0th, 1st adn 2nd mass moments of the droplet or
 raindroplet distributions after microphysics or motion */
template <typename MonitorViewsType>
struct MonitorMassMoments {
  MonitorViewsType microphysics_monitor;  // monitoring mass moments during microphysics
  MonitorViewsType motion_monitor;        // monitoring mass moments during motion

  /**
   * @brief Reset monitors for mass moments from both motion and microphysics.
   */
  void reset_monitor() const {
    microphysics_monitor.reset_views();
    motion_monitor.reset_views();
  }

  /**
   * @brief Placeholder function to obey SDMMonitor concept does nothing.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param totmass_condensed Mass condensed in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_condensation(const TeamMember& team_member, const double totmass_condensed) const {}

  /**
   * @brief Monitor 0th, 1st and 2nd moments of the droplet mass distribution
   *
   * calls calculate_massmoments to monitor the moments of the droplet mass
   * distribution during SDM microphysics
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember& team_member, const viewd_constsupers supers) const {
    Kokkos::single(Kokkos::PerTeam(team_member), [=, this]() {
      microphysics_monitor.calculate_massmoments(team_member, supers);
    });
  }

  /**
   * @brief Monitor 0th, 1st and 2nd moments of the droplet mass distribution
   *
   * calls calculate_massmoments to monitor the moments of the droplet mass
   * distribution during SDM motion.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one motion timestep
   */
  void monitor_motion(const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for(
        "monitor_motion", TeamPolicy(ngbxs, Kokkos::AUTO()),
        KOKKOS_CLASS_LAMBDA(const TeamMember& team_member) {
          Kokkos::single(Kokkos::PerTeam(team_member), [=, this]() {
            const auto ii = team_member.league_rank();
            const auto supers(d_gbxs(ii).supersingbx.readonly());
            motion_monitor.calculate_massmoments(team_member, supers);
          });
        });
  }

  /**
   * @brief Constructor for MonitorMassMoments
   *
   * @param ngbxs Number of gridboxes in domain.
   */
  explicit MonitorMassMoments(const size_t ngbxs)
      : microphysics_monitor(ngbxs), motion_monitor(ngbxs) {
    reset_monitor();
  }
};

#endif  //  LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_HPP_
