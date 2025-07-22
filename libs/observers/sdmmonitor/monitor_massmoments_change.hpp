/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_massmoments_change.hpp
 * Project: sdmmonitor
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct to create observer which outputs the change in mass moments due to motion and microphysics
 * (seperately) in each gridbox over a constant interval (output at the start of each timestep).
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_CHANGE_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_CHANGE_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <cstdint>
#include <memory>

#include "../../kokkosaliases.hpp"
#include "observers/massmoments_observer.hpp"
#include "zarr/buffer.hpp"

namespace KCS = KokkosCleoSettings;

struct MonitorMassMomentsChangeViews {
  Buffer<uint64_t>::mirrorviewd_buffer
      d_mom0_prev;  // view on device for storing previous 0th mass moment
  Buffer<float>::mirrorviewd_buffer
      d_mom1_prev;  // view on device for storing previous 1st mass moment
  Buffer<float>::mirrorviewd_buffer
      d_mom2_prev;  // view on device for storing previous 2nd mass moment

  Buffer<uint64_t>::mirrorviewd_buffer
      d_delta_mom0;  // view on device for monitoring change in 0th mass moment
  Buffer<float>::mirrorviewd_buffer
      d_delta_mom1;  // view on device for monitoring change in 1st mass moment
  Buffer<float>::mirrorviewd_buffer
      d_delta_mom2;  // view on device for monitoring change in 2nd mass moment

  /**
   * @brief Parallel loop to fill device views for change in mass moments with zero value
   */
  void reset_views() const {
    Kokkos::parallel_for(
        "reset_views", Kokkos::RangePolicy(0, d_mom0.extent(0)),
        KOKKOS_CLASS_LAMBDA(const size_t jj) {
          d_delta_mom0(jj) = 0;
          d_delta_mom1(jj) = 0.0;
          d_delta_mom2(jj) = 0.0;
        });
  }

  /**
   * @brief Write the change in 0th, 1st and 2nd moments of the droplet mass distribution to views.
   *
   * Calculates the current mass moments and writes the change in their values since they were last
   * calculated to the mass moment deltas (d_delta_mom0, d_delta_mom1 and d_delta_mom2).
   *
   * _Note:_ possible conversion of mass moments at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox
   */
  KOKKOS_FUNCTION
  void fetch_delta_massmoments(const TeamMember& team_member,
                               const viewd_constsupers supers) const {
    const auto ii = team_member.league_rank();
    auto mom0_now = 0;
    auto mom1_now = 0.0;
    auto mom2_now = 0.0;
    calculate_massmoments(team_member, supers, mom0_now, mom1_now, mom2_now);

    /* accumulate change in mass moments */
    d_delta_mom0 += mom0_now - d_mom0_prev(ii);
    d_delta_mom1 += mom1_now - d_mom1_prev(ii);
    d_delta_mom2 += mom2_now - d_mom2_prev(ii);

    /* store current mass moments as previous one for next accumulation */
    d_mom0_prev(ii) = mom0_now;
    d_mom1_prev(ii) = mom1_now;
    d_mom2_prev(ii) = mom2_now;
  }

  explicit MonitorMassMomentsChangeViews(const size_t ngbxs)
      : d_mom0_prev("d_monitor_mom0_prev", ngbxs),
        d_mom1_prev("d_monitor_mom1_prev", ngbxs),
        d_mom2_prev("d_monitor_mom2_prev", ngbxs),
        d_delta_mom0("d_monitor_delta_mom0", ngbxs),
        d_delta_mom1("d_monitor_delta_mom1", ngbxs),
        d_delta_mom2("d_monitor_delta_mom2", ngbxs) {
    reset_views();
  }
};

struct MonitorRainMassMomentsChangeViews {
  Buffer<uint64_t>::mirrorviewd_buffer
      d_mom0_prev;  // view on device for storing previous 0th mass moment
  Buffer<float>::mirrorviewd_buffer
      d_mom1_prev;  // view on device for storing previous 1st mass moment
  Buffer<float>::mirrorviewd_buffer
      d_mom2_prev;  // view on device for storing previous 2nd mass moment

  Buffer<uint64_t>::mirrorviewd_buffer
      d_delta_mom0;  // view on device for monitoring change in 0th mass moment
  Buffer<float>::mirrorviewd_buffer
      d_delta_mom1;  // view on device for monitoring change in 1st mass moment
  Buffer<float>::mirrorviewd_buffer
      d_delta_mom2;  // view on device for monitoring change in 2nd mass moment

  /**
   * @brief Parallel loop to fill device views with zero value
   */
  void reset_views() const {
    Kokkos::parallel_for(
        "reset_views", Kokkos::RangePolicy(0, d_mom0.extent(0)),
        KOKKOS_CLASS_LAMBDA(const size_t jj) {
          d_delta_mom0(jj) = 0;
          d_delta_mom1(jj) = 0.0;
          d_delta_mom2(jj) = 0.0;
        });
  }

  /**
   * @brief Write the 0th, 1st and 2nd moments of the raindroplet mass distribution to data views.
   *
   * Calculates the current mass moments of the raindrop distribution and overwrites the current
   * values for the raindrop mass moments (d_mom0, d_mom1 and d_mom2) stored since the data
   * views were last reset.
   *
   * _Note:_ possible conversion of raindrop mass moments at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox
   */
  KOKKOS_FUNCTION
  void fetch_delta_massmoments(const TeamMember& team_member,
                               const viewd_constsupers supers) const {
    const auto ii = team_member.league_rank();
    auto mom0_now = 0;
    auto mom1_now = 0.0;
    auto mom2_now = 0.0;
    calculate_rainmassmoments(team_member, supers, mom0_now, mom1_now, mom2_now);

    /* accumulate change in mass moments */
    d_delta_mom0 += mom0_now - d_mom0_prev(ii);
    d_delta_mom1 += mom1_now - d_mom1_prev(ii);
    d_delta_mom2 += mom2_now - d_mom2_prev(ii);

    /* store current mass moments as previous one for next accumulation */
    d_mom0_prev(ii) = mom0_now;
    d_mom1_prev(ii) = mom1_now;
    d_mom2_prev(ii) = mom2_now;
  }

  explicit MonitorRainMassMomentsChangeViews(const size_t ngbxs)
      : d_mom0_prev("d_monitor_rainmom0_prev", ngbxs),
        d_mom1_prev("d_monitor_rainmom1_prev", ngbxs),
        d_mom2_prev("d_monitor_rainmom2_prev", ngbxs),
        d_delta_mom0("d_monitor_rain_delta_mom0", ngbxs),
        d_delta_mom1("d_monitor_rain_delta_mom1", ngbxs),
        d_delta_mom2("d_monitor_rain_delta_mom2", ngbxs) {
    reset_views();
  }
};

/* struct satisfies SDMMonitor concept in order to make observer for monitoring mass moments
 * according to the templated MonitorViewsType e.g. 0th, 1st adn 2nd mass moments of the droplet or
 raindroplet distribution after microphysics and/or motion */
template <typename MonitorViewsType>
struct MonitorMassMomentsChange {
  MonitorViewsType microphysics_moms;  // mass moments monitored during microphysics
  MonitorViewsType motion_moms;        // mass moments monitored during motion

  /**
   * @brief Reset monitors for mass moments from both motion and microphysics.
   */
  void reset_monitor() const {
    microphysics_moms.reset_views();
    motion_moms.reset_views();
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
   * calls fetch_delta_massmoments to monitor the moments of the droplet mass
   * distribution during SDM microphysics
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember& team_member, const viewd_constsupers supers) const {
    microphysics_moms.fetch_delta_massmoments(team_member, supers);
  }

  /**
   * @brief Monitor 0th, 1st and 2nd moments of the droplet mass distribution
   *
   * calls fetch_delta_massmoments to monitor the moments of the droplet mass
   * distribution during SDM motion.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   * @param domainsupers The view of all super-droplets (in bounds of domain).
   */
  void monitor_motion(const viewd_constgbx d_gbxs, const subviewd_constsupers domainsupers) const {
    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for(
        "monitor_motion", TeamPolicy(ngbxs, KCS::team_size),
        KOKKOS_CLASS_LAMBDA(const TeamMember& team_member) {
          const auto ii = team_member.league_rank();
          const auto supers = d_gbxs(ii).supersingbx.readonly(domainsupers);
          motion_moms.fetch_delta_massmoments(team_member, supers);
        });
  }

  /**
   * @brief Constructor for MonitorMassMomentsChange
   *
   * @param ngbxs Number of gridboxes in domain.
   */
  explicit MonitorMassMomentsChange(const size_t ngbxs)
      : microphysics_moms(ngbxs), motion_moms(ngbxs) {
    reset_monitor();
  }
};

#endif  // LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_CHANGE_HPP_
