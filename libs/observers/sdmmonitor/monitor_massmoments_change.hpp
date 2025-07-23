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
#include "gridboxes/gridboxmaps.hpp"
#include "observers/massmoments_observer.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr/buffer.hpp"

namespace KCS = KokkosCleoSettings;

struct MonitorMassMomentsChangeViews {
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
        "reset_views", Kokkos::RangePolicy(0, d_delta_mom0.extent(0)),
        KOKKOS_CLASS_LAMBDA(const size_t jj) {
          d_delta_mom0(jj) = 0;
          d_delta_mom1(jj) = 0.0;
          d_delta_mom2(jj) = 0.0;
        });
  }

  /**
   * @brief Before timestepping write the 0th, 1st and 2nd moments of the raindroplet mass
   * distribution to "prev" views.
   *
   * Calculates the current mass moments of the raindrop distribution and stores them in the
   * "prev" views (d_mom0_prev, d_mom1_prev and d_mom2_prev), so the change in the moments can be
   * calculated by fetch_delta_massmoments during the first timestep (and onwards).
   *
   * _Note:_ possible conversion of raindrop mass moments at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox
   * @param d_mom0_prev View on device of previous 0th mass moment
   * @param d_mom1_prev View on device of previous 1th mass moment
   * @param d_mom2_prev View on device of previous 2th mass moment
   */
  KOKKOS_FUNCTION
  void before_timestepping(const TeamMember& team_member, const viewd_constsupers supers,
                           Buffer<uint64_t>::mirrorviewd_buffer d_mom0_prev,
                           Buffer<float>::mirrorviewd_buffer d_mom1_prev,
                           Buffer<float>::mirrorviewd_buffer d_mom2_prev) const {
    const auto ii = team_member.league_rank();
    calculate_massmoments(team_member, supers, d_mom0_prev(ii), d_mom1_prev(ii), d_mom2_prev(ii));
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
   * @param d_mom0_prev View on device of previous 0th mass moment
   * @param d_mom1_prev View on device of previous 1th mass moment
   * @param d_mom2_prev View on device of previous 2th mass moment
   */
  KOKKOS_FUNCTION
  void fetch_delta_massmoments(const TeamMember& team_member, const viewd_constsupers supers,
                               Buffer<uint64_t>::mirrorviewd_buffer d_mom0_prev,
                               Buffer<float>::mirrorviewd_buffer d_mom1_prev,
                               Buffer<float>::mirrorviewd_buffer d_mom2_prev) const {
    const auto ii = team_member.league_rank();

    uint64_t mom0_now = 0;
    float mom1_now = 0.0;
    float mom2_now = 0.0;
    calculate_massmoments(team_member, supers, mom0_now, mom1_now, mom2_now);

    /* accumulate change in mass moments */
    d_delta_mom0(ii) += mom0_now - d_mom0_prev(ii);
    d_delta_mom1(ii) += mom1_now - d_mom1_prev(ii);
    d_delta_mom2(ii) += mom2_now - d_mom2_prev(ii);

    /* store current mass moments as previous one for next accumulation */
    d_mom0_prev(ii) = mom0_now;
    d_mom1_prev(ii) = mom1_now;
    d_mom2_prev(ii) = mom2_now;
  }

  explicit MonitorMassMomentsChangeViews(const size_t ngbxs)
      : d_delta_mom0("d_monitor_delta_mom0", ngbxs),
        d_delta_mom1("d_monitor_delta_mom1", ngbxs),
        d_delta_mom2("d_monitor_delta_mom2", ngbxs) {
    reset_views();
  }
};

struct MonitorRainMassMomentsChangeViews {
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
        "reset_views", Kokkos::RangePolicy(0, d_delta_mom0.extent(0)),
        KOKKOS_CLASS_LAMBDA(const size_t jj) {
          d_delta_mom0(jj) = 0;
          d_delta_mom1(jj) = 0.0;
          d_delta_mom2(jj) = 0.0;
        });
  }

  /**
   * @brief Before timestepping write the 0th, 1st and 2nd moments of the raindroplet mass
   * distribution to "prev" views.
   *
   * Calculates the current mass moments of the raindrop distribution and stores them in the
   * "prev" views (d_mom0_prev, d_mom1_prev and d_mom2_prev), so the change in the moments can be
   * calculated by fetch_delta_massmoments during the first timestep (and onwards).
   *
   * _Note:_ possible conversion of raindrop mass moments at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox
   * @param d_mom0_prev View on device of previous 0th mass moment
   * @param d_mom1_prev View on device of previous 1th mass moment
   * @param d_mom2_prev View on device of previous 2th mass moment
   */
  KOKKOS_FUNCTION
  void before_timestepping(const TeamMember& team_member, const viewd_constsupers supers,
                           Buffer<uint64_t>::mirrorviewd_buffer d_mom0_prev,
                           Buffer<float>::mirrorviewd_buffer d_mom1_prev,
                           Buffer<float>::mirrorviewd_buffer d_mom2_prev) const {
    const auto ii = team_member.league_rank();
    calculate_rainmassmoments(team_member, supers, d_mom0_prev(ii), d_mom1_prev(ii),
                              d_mom2_prev(ii));
  }

  /**
   * @brief Write the change in the 0th, 1st and 2nd moments of the raindroplet mass distribution
   * to data views.
   *
   * Calculates the current mass moments of the raindrop distribution and writes the change in their
   * values since they were last calculated to the mass moment deltas
   * (d_delta_mom0, d_delta_mom1 and d_delta_mom2).
   *
   * _Note:_ possible conversion of raindrop mass moments at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox
   * @param d_mom0_prev View on device of previous 0th mass moment
   * @param d_mom1_prev View on device of previous 1th mass moment
   * @param d_mom2_prev View on device of previous 2th mass moment
   */
  KOKKOS_FUNCTION
  void fetch_delta_massmoments(const TeamMember& team_member, const viewd_constsupers supers,
                               Buffer<uint64_t>::mirrorviewd_buffer d_mom0_prev,
                               Buffer<float>::mirrorviewd_buffer d_mom1_prev,
                               Buffer<float>::mirrorviewd_buffer d_mom2_prev) const {
    const auto ii = team_member.league_rank();

    uint64_t mom0_now = 0;
    float mom1_now = 0.0;
    float mom2_now = 0.0;
    calculate_rainmassmoments(team_member, supers, mom0_now, mom1_now, mom2_now);

    /* accumulate change in mass moments */
    d_delta_mom0(ii) += mom0_now - d_mom0_prev(ii);
    d_delta_mom1(ii) += mom1_now - d_mom1_prev(ii);
    d_delta_mom2(ii) += mom2_now - d_mom2_prev(ii);

    /* store current mass moments as previous one for next accumulation */
    d_mom0_prev(ii) = mom0_now;
    d_mom1_prev(ii) = mom1_now;
    d_mom2_prev(ii) = mom2_now;
  }

  explicit MonitorRainMassMomentsChangeViews(const size_t ngbxs)
      : d_delta_mom0("d_monitor_rain_delta_mom0", ngbxs),
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
  Buffer<uint64_t>::mirrorviewd_buffer
      d_mom0_prev;  // view on device for storing previous 0th mass moment
  Buffer<float>::mirrorviewd_buffer
      d_mom1_prev;  // view on device for storing previous 1st mass moment
  Buffer<float>::mirrorviewd_buffer
      d_mom2_prev;  // view on device for storing previous 2nd mass moment

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
   * @param d_gbxs The view of gridboxes in device memory.
   */
  KOKKOS_FUNCTION
  void before_timestepping(const TeamMember& team_member,
                           const subviewd_constsupers d_supers) const {
    motion_moms.before_timestepping(
        team_member, d_supers, d_mom0_prev, d_mom1_prev,
        d_mom2_prev);  // same outcome as microphysics_moms.before_timestepping(...);
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
    microphysics_moms.fetch_delta_massmoments(team_member, supers, d_mom0_prev, d_mom1_prev,
                                              d_mom2_prev);
  }

  /**
   * @brief Monitor 0th, 1st and 2nd moments of the droplet mass distribution
   *
   * calls fetch_delta_massmoments to monitor the moments of the droplet mass
   * distribution during SDM motion
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_motion(const TeamMember& team_member, const viewd_constsupers supers) const {
    motion_moms.fetch_delta_massmoments(team_member, supers, d_mom0_prev, d_mom1_prev, d_mom2_prev);
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
        "monitor_motion_massmoments", TeamPolicy(ngbxs, KCS::team_size),
        KOKKOS_CLASS_LAMBDA(const TeamMember& team_member) {
          const auto ii = team_member.league_rank();
          const auto supers = d_gbxs(ii).supersingbx.readonly(domainsupers);
          monitor_motion(team_member, supers);
        });
  }

  /**
   * @brief Placeholder function to obey SDMMonitor concept does nothing.
   *
   * @param gbxindex gridbox whose bottom boundary is to be evaluated.
   * @param gbxmaps The Gridbox Maps.
   * @param state The State of the volume containing the super-droplets (gridbox matching gbxindex).
   * @param drop The super-droplet to evaluate.
   */
  KOKKOS_FUNCTION
  void monitor_precipitation(const TeamMember& team_member, const unsigned int gbxindex,
                             const GridboxMaps auto& gbxmaps, Superdrop& drop) const {}

  /**
   * @brief Constructor for MonitorMassMomentsChange
   *
   * @param ngbxs Number of gridboxes in domain.
   */
  explicit MonitorMassMomentsChange(const size_t ngbxs)
      : microphysics_moms(ngbxs),
        motion_moms(ngbxs),
        d_mom0_prev("d_monitor_mom0_prev", ngbxs),
        d_mom1_prev("d_monitor_mom1_prev", ngbxs),
        d_mom2_prev("d_monitor_mom2_prev", ngbxs) {
    reset_monitor();
  }
};

#endif  // LIBS_OBSERVERS_SDMMONITOR_MONITOR_MASSMOMENTS_CHANGE_HPP_
