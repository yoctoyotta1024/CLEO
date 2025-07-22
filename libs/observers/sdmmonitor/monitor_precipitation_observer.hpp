/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_precipitation_observer.hpp
 * Project: sdmmonitor
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct to create observer which outputs accumulated precipitation over a constant timestep, i.e.
 * the mean rate of precipitation over a timestep, by monitoring superdroplet motion through
 * the bottom boundary of each gridbox,
 * i.e. output = downward mass flux of water / water density * timestep
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_MONITOR_PRECIPITATION_OBSERVER_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_MONITOR_PRECIPITATION_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>

#include "../../cleoconstants.hpp"
#include "../../kokkosaliases.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/observers.hpp"
#include "observers/sdmmonitor/do_sdmmonitor_obs.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

namespace dlc = dimless_constants;

/* struct satisfies SDMMonitor concept for use in do_sdmmonitor_obs to make observer */
struct MonitorPrecipitation {
  using datatype = double;
  Buffer<datatype>::mirrorviewd_buffer d_data;  // view on device copied to host by DoSDMMonitorObs

  /**
   * @brief Parallel loop to fill d_data with zero value.
   */
  void reset_monitor() const;

  /**
   * @brief Placeholder function to obey SDMMonitor concept does nothing.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   */
  KOKKOS_FUNCTION
  void before_timestepping(const TeamMember& team_member,
                           const subviewd_constsupers d_supers) const {}

  /**
   * @brief Placeholder function to obey SDMMonitor concept does nothing.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param totmass_condensed Mass condensed in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_condensation(const TeamMember& team_member, const double totmass_condensed) const {}

  /**
   * @brief Placeholder function to obey SDMMonitor concept does nothing.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param supers (sub)View of all the superdrops in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember& team_member, const viewd_constsupers supers) const {}

  /**
   * @brief Placeholder function to obey SDMMonitor concept does nothing.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   * @param domainsupers The view of superdroplets within the domain in device memory.
   */
  void monitor_motion(const viewd_constgbx d_gbxs, const subviewd_constsupers domainsupers) const {}

  /**
   * @brief calculate accumulated precipitation over a constant timestep, i.e.
   * the mean rate of precipitation over a timestep, as the droplet motion through
   * the bottom boundary of each gridbox,
   * i.e. output = downward mass flux of water / water density * timestep
   *
   * @param gbxindex gridbox whose bottom boundary is to be evaluated.
   * @param gbxmaps The Gridbox Maps.
   * @param state The State of the volume containing the super-droplets (gridbox matching gbxindex).
   * @param drop The super-droplet to evaluate.
   */
  KOKKOS_FUNCTION
  void monitor_precipitation(const TeamMember& team_member, const unsigned int gbxindex,
                             const GridboxMaps auto& gbxmaps, Superdrop& drop) const {
    const auto lowerlim = gbxmaps.coord3bounds(gbxindex).first;
    if (drop.get_coord3() < lowerlim) {
      const auto ii = team_member.league_rank();
      const auto gbxarea = gbxmaps.get_gbxarea(gbxindex);
      const auto m = drop.condensate_mass() * drop.get_xi() / dlc::Rho_l / gbxarea;
      Kokkos::atomic_add(&d_data(ii), m);
    }
  }

  /**
   * @brief Constructor for MonitorPrecipitation
   *
   * @param ngbxs Number of gridboxes in domain.
   */
  explicit MonitorPrecipitation(const size_t ngbxs) : d_data("precip", ngbxs) { reset_monitor(); }
};

/**
 * @brief Constructs an observer which writes data monitoring precipitaion (i.e. the downward mass
 * flux of water / water density * timestep) to an array with a constant observation timestep
 * "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Constructed type satisfying observer concept.
 */
template <typename Dataset, typename Store>
inline Observer auto MonitorPrecipitationObserver(const unsigned int interval, Dataset& dataset,
                                                  Store& store, const size_t maxchunk,
                                                  const size_t ngbxs) {
  using Mo = MonitorPrecipitation;
  const auto name = std::string_view("precip");
  const auto units = std::string_view("m");
  constexpr auto scale_factor = dlc::R0 * dlc::R0 * dlc::R0 / dlc::COORD0 / dlc::COORD0;
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr_ptr = std::make_shared<XarrayZarrArray<Store, Mo::datatype>>(
      dataset.template create_array<Mo::datatype>(name, units, scale_factor, chunkshape, dimnames));

  const auto do_obs =
      DoSDMMonitorObs<Dataset, Store, Mo, Mo::datatype>(dataset, store, xzarr_ptr, Mo(ngbxs));
  return ConstTstepObserver(interval, do_obs);
}

#endif  // LIBS_OBSERVERS_SDMMONITOR_MONITOR_PRECIPITATION_OBSERVER_HPP_
