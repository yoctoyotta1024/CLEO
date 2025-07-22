/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_condensation_observer.hpp
 * Project: sdmmonitor
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct to create observer which outputs mass change due to condensation monitored from SDM
 * microphysical process in each gridbox a constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_MONITOR_CONDENSATION_OBSERVER_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_MONITOR_CONDENSATION_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>

#include "../../kokkosaliases.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/observers.hpp"
#include "observers/sdmmonitor/do_sdmmonitor_obs.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

/* struct satisfies SDMMonitor concept for use in do_sdmmonitor_obs to make observer */
struct MonitorCondensation {
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
   * @brief Monitor mass of liquid change due to condensation / evaporation
   *
   * Add totmass_condensed to current value for mass condensed since d_data was last reset.
   *
   * _Note:_ possible conversion of mass condensed at one timestep from double precision
   * (8 bytes double) to single precision (4 bytes float) in output depending on datatype alias.
   *
   * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
   * @param totmass_condensed Mass condensed in one gridbox during one microphysical timestep
   */
  KOKKOS_FUNCTION
  void monitor_condensation(const TeamMember& team_member, const double totmass_condensed) const;

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
   * @brief Constructor for MonitorCondensation
   *
   * @param ngbxs Number of gridboxes in domain.
   */
  explicit MonitorCondensation(const size_t ngbxs) : d_data("massdelta_cond", ngbxs) {
    reset_monitor();
  }
};

/**
 * @brief Constructs an observer which writes data monitoring condensation microphysics to an
 * array with a constant observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Constructed type satisfying observer concept.
 */
template <typename Dataset, typename Store>
inline Observer auto MonitorCondensationObserver(const unsigned int interval, Dataset& dataset,
                                                 Store& store, const size_t maxchunk,
                                                 const size_t ngbxs) {
  using Mo = MonitorCondensation;
  const auto name = std::string_view("massdelta_cond");
  const auto units = std::string_view("g");
  constexpr auto scale_factor = dlc::MASS0grams;
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr_ptr = std::make_shared<XarrayZarrArray<Store, Mo::datatype>>(
      dataset.template create_array<Mo::datatype>(name, units, scale_factor, chunkshape, dimnames));

  const auto do_obs =
      DoSDMMonitorObs<Dataset, Store, Mo, Mo::datatype>(dataset, store, xzarr_ptr, Mo(ngbxs));
  return ConstTstepObserver(interval, do_obs);
}

#endif  // LIBS_OBSERVERS_SDMMONITOR_MONITOR_CONDENSATION_OBSERVER_HPP_
