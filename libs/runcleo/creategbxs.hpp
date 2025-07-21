/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: creategbxs.hpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * classes and templated functionality for creating initialised view of Gridboxes from some
 * initial conditions
 */

#ifndef LIBS_RUNCLEO_CREATEGBXS_HPP_
#define LIBS_RUNCLEO_CREATEGBXS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../kokkosaliases.hpp"
#include "gridboxes/findrefs.hpp"
#include "gridboxes/gbxindex.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/supersindomain.hpp"
#include "gridboxes/supersingbx.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

namespace KCS = KokkosCleoSettings;

/**
 * @brief A generator for initial Gridboxes.
 *
 * This class provides functionality to generate Gridbox instances based on
 * the gridbox maps and on the initial conditions stored in this struct's vectors.
 */
class GenGridbox {
 private:
  std::shared_ptr<Gbxindex::Gen> GbxindexGen;
  /**< Pointer to gridbox index generator, Gbxindex::Gen object */
  std::vector<double> presss; /**< Vector of pressures for each gridbox */
  std::vector<double> temps;  /**< Vector of temperatures for each gridbox */
  std::vector<double> qvaps;  /**< Vector of vapor mass mixing ratio for each gridbox */
  std::vector<double> qconds; /**< Vector of condensed water mass mixing ratio for each gridbox */
  std::vector<std::pair<double, double>> wvels;
  /**< Vector of vertical (coord3) wind velocities for each gridbox */
  std::vector<std::pair<double, double>> uvels;
  /**< Vector of eastward (coord1) wind velocities for each gridbox */
  std::vector<std::pair<double, double>> vvels;
  /**< Vector of northward (coord2) wind velocities for each gridbox */

  /**
   * @brief Get the state of a specified Gridbox from the initial conditions.
   *
   * This function returns the state of the Gridbox at the ii'th index in the initial conditions
   * given by the GenGridbox struct.
   *
   * @param ii The index of the Gridbox in the initial conditions data.
   * @param volume The volume of the Gridbox.
   * @return The State of the Gridbox.
   */
  State state_at(const unsigned int ii, const double volume) const;

 public:
  /**
   * @brief Constructs a GenGridbox object.
   *
   * Constructs a GenGridbox object based on the provided initial conditions in 'GbxInitConds'.
   *
   * @tparam GbxInitConds Type of the Gridboxes' initial conditions.
   * @param gbxic The initial conditions for the Gridboxes.
   */
  template <typename GbxInitConds>
  explicit GenGridbox(const GbxInitConds &gbxic)
      : GbxindexGen(std::make_shared<Gbxindex::Gen>()),
        presss(gbxic.press()),
        temps(gbxic.temp()),
        qvaps(gbxic.qvap()),
        qconds(gbxic.qcond()),
        wvels(gbxic.wvel()),
        uvels(gbxic.uvel()),
        vvels(gbxic.vvel()) {}

  /**
   * @brief Serial version of operator to generate a Gridbox from the data at the ii'th position
   * of the initial conditions data.
   *
   * This function generates a Gridbox corresponding to the ii'th position in the
   * initial conditions data and using the gridbox maps and view of (all) superdroplets.
   *
   * @tparam GbxMaps Type of the Gridbox Maps.
   * @param ii The index of the Gridbox.
   * @param gbxmaps The Gridbox Maps.
   * @param domainsupers The view of all super-droplets (in bounds of domain).
   * @return The generated Gridbox.
   */
  template <GridboxMaps GbxMaps>
  Gridbox operator()(const unsigned int ii, const GbxMaps &gbxmaps,
                     const subviewd_constsupers domainsupers) const {
    const auto gbxindex = GbxindexGen->next(ii);
    const auto volume = gbxmaps.get_gbxvolume(gbxindex.value);
    const auto state = state_at(ii, volume);

    return Gridbox(gbxindex, state, domainsupers);
  }

  /**
   * @brief Parallel-safe version of operator to generate a Gridbox from the data at the
   * ii'th position of the initial conditions data.
   *
   * This function generates a gridbox at the specified index using gridbox maps,
   * total superdroplets, and host team member.
   * Given a Kokkos team thread ('team_member'), this function generates a Gridbox corresponding to
   * the ii'th position in the initial conditions data and using the gridbox maps and view of
   * (all) superdroplets on device and host.
   *
   * @tparam GbxMaps Type of the Gridbox Maps.
   * @param team_member The host team member reference.
   * @param ii The index of the Gridbox.
   * @param gbxmaps The Gridbox Maps.
   * @param h_domainsupers Host mirror of view of all super-droplets (in bounds of domain).
   * @return The generated Gridbox.
   */
  template <GridboxMaps GbxMaps>
  Gridbox operator()(const HostTeamMember &team_member, const unsigned int ii,
                     const GbxMaps &gbxmaps,
                     const subviewd_constsupers::HostMirror h_domainsupers) const {
    const auto gbxindex = GbxindexGen->next(ii);
    const auto volume = gbxmaps.get_gbxvolume(gbxindex.value);
    const auto state = state_at(ii, volume);
    const auto refs = find_refs(team_member, h_domainsupers, gbxindex.value);

    return Gridbox(gbxindex, state, refs);
  }
};

/**
 * @brief Create Gridboxes from initial conditions
 *
 * This function creates Gridboxes based on the provided gridbox maps and initial conditions,
 * and given super-droplets.
 *
 * @tparam GbxMaps Type representing Gridbox Maps.
 * @tparam GbxInitConds Type representing Gridbox initial conditions.
 *
 * @param gbxmaps The Gridbox Maps.
 * @param gbxic The Gridbox initial conditions.
 * @param allsupers The view of all super-droplets (in bounds of domain).
 *
 * @return The view of initialised Gridboxes.
 */
template <GridboxMaps GbxMaps, typename GbxInitConds>
dualview_gbx create_gbxs(const GbxMaps &gbxmaps, const GbxInitConds &gbxic,
                         const SupersInDomain &allsupers);

/**
 * @brief Initialise host view of Gridboxes.
 *
 * This function initialises Gridboxes in host memory using data from Gridbox Maps,
 * a Gridbox generator, and a view of super-droplets.
 *
 * The equivalent serial version of Kokkos::parallel_for([...]) loop is:
 * @code
 * for (size_t ii(0); ii < ngbxs; ++ii)
 * {
 *  h_gbxs(ii) = gen(ii, gbxmaps, domainsupers);
 * }
 * @endcode
 *
 * @tparam GbxMaps Type representing Gridbox Maps.
 *
 * @param gbxmaps The Gridbox Maps.
 * @param gen The Gridbox generator.
 * @param domainsupers The view of all super-droplets (in bounds of domain).
 * @param h_gbxs The view of Gridboxes on the host.
 */
template <GridboxMaps GbxMaps>
inline void initialise_gbxs_on_host(const GbxMaps &gbxmaps, const GenGridbox &gen,
                                    const subviewd_constsupers domainsupers,
                                    const viewh_gbx h_gbxs);

/**
 * @brief Initialise a view of Gridboxes.
 *
 * This function initialises a view of gridboxes in device memory using data
 * from an instance of GbxInitConds for each Gridbox's index, initial State and
 * view of super-droplets.
 *
 * @tparam GbxMaps Type representing Gridbox Maps.
 * @tparam GbxInitConds Type representing Gridbox initial conditions.
 *
 * @param gbxmaps The Gridbox Maps.
 * @param gbxic The initial conditions for the Gridboxes.
 * @param domainsupers The view of all super-droplets (in bounds of domain).
 *
 * @return The initialised view of Gridboxes.
 */
template <GridboxMaps GbxMaps, typename GbxInitConds>
inline dualview_gbx initialise_gbxs(const GbxMaps &gbxmaps, const GbxInitConds &gbxic,
                                    const subviewd_constsupers domainsupers);

/**
 * @brief Check if gridbox initialisation is complete.
 *
 * This function checks if the number of created gridboxes is consistent with
 * the number of gridboxes from gridbox maps and if each gridbox has correct
 * references to superdroplets.
 *
 * @param ngbxs_from_maps The number of gridboxes from gridbox maps.
 * @param gbxs The dualview containing gridboxes.
 * @param totsupers The view containing superdroplets on device.
 */
void is_gbxinit_complete(const size_t ngbxs_from_maps, dualview_gbx gbxs,
                         const viewd_constsupers totsupers);

/**
 * @brief Print some information about initial Gridboxes.
 *
 * This function prints information about each Gridbox, including its index,
 * volume, and number of super-droplets.
 *
 * @param h_gbxs The host view of Gridboxes.
 */
void print_gbxs(const viewh_constgbx gbxs);

/**
 * @brief Create Gridboxes from initial conditions
 *
 * This function creates Gridboxes based on the provided gridbox maps and initial conditions,
 * and given super-droplets.
 *
 * Kokkos::Profiling are null pointers unless a Kokkos profiler library has been
 * exported to "KOKKOS_TOOLS_LIBS" prior to runtime so the lib gets dynamically loaded.
 *
 * @tparam GbxMaps Type representing Gridbox Maps.
 * @tparam GbxInitConds Type representing Gridbox initial conditions.
 *
 * @param gbxmaps The Gridbox Maps.
 * @param gbxic The Gridbox initial conditions.
 * @param allsupers The view of all super-droplets (in bounds of domain).
 *
 * @return The view of initialised Gridboxes.
 */
template <GridboxMaps GbxMaps, typename GbxInitConds>
dualview_gbx create_gbxs(const GbxMaps &gbxmaps, const GbxInitConds &gbxic,
                         const SupersInDomain &allsupers) {
  Kokkos::Profiling::ScopedRegion region("init_gbxs");

  std::cout << "\n--- create gridboxes ---\ninitialising\n";
  const auto domainsupers = allsupers.domain_supers_readonly();
  const auto gbxs = initialise_gbxs(gbxmaps, gbxic, domainsupers);

  std::cout << "checking initialisation\n";
  is_gbxinit_complete(gbxmaps.get_local_ngridboxes_hostcopy(), gbxs,
                      allsupers.get_totsupers_readonly());

  // // Print information about the created superdrops
  // print_gbxs(gbxs.view_host());

  std::cout << "--- create gridboxes: success ---\n";

  return gbxs;
}

/**
 * @brief Initialise a view of Gridboxes.
 *
 * This function initialises a view of gridboxes in device memory using data
 * from an instance of GbxInitConds for each Gridbox's index, initial State and
 * view of super-droplets.
 *
 * @tparam GbxMaps Type representing Gridbox Maps.
 * @tparam GbxInitConds Type representing Gridbox initial conditions.
 *
 * @param gbxmaps The Gridbox Maps.
 * @param gbxic The initial conditions for the Gridboxes.
 * @param domainsupers The view of all super-droplets (in bounds of domain).
 *
 * @return The initialised view of Gridboxes.
 */
template <GridboxMaps GbxMaps, typename GbxInitConds>
inline dualview_gbx initialise_gbxs(const GbxMaps &gbxmaps, const GbxInitConds &gbxic,
                                    const subviewd_constsupers domainsupers) {
  // create dualview for gridboxes on device and host memory
  dualview_gbx gbxs("gbxs", gbxic.get_ngbxs());

  // initialise gridboxes on host
  const GenGridbox gen(gbxic);
  gbxs.sync_host();
  initialise_gbxs_on_host(gbxmaps, gen, domainsupers, gbxs.view_host());
  gbxs.modify_host();

  // update device gridbox view to match host's gridbox view
  gbxs.sync_device();

  return gbxs;
}

/**
 * @brief Initialise host view of Gridboxes.
 *
 * This function initialises Gridboxes in host memory using data from Gridbox Maps,
 * a Gridbox generator, and a view of super-droplets.
 *
 * The equivalent serial version of Kokkos::parallel_for([...]) loop is:
 * @code
 * for (size_t ii(0); ii < ngbxs; ++ii)
 * {
 *  h_gbxs(ii) = gen(ii, gbxmaps, domainsupers);
 * }
 * @endcode
 *
 * @tparam GbxMaps Type representing Gridbox Maps.
 *
 * @param gbxmaps The Gridbox Maps.
 * @param gen The Gridbox generator.
 * @param domainsupers The view of all super-droplets (both in and out of bounds of domain).
 * @param h_gbxs The view of Gridboxes on the host.
 */
template <GridboxMaps GbxMaps>
inline void initialise_gbxs_on_host(const GbxMaps &gbxmaps, const GenGridbox &gen,
                                    const subviewd_constsupers domainsupers,
                                    const viewh_gbx h_gbxs) {
  const size_t ngbxs = h_gbxs.extent(0);

  auto h_domainsupers = Kokkos::create_mirror_view(
      domainsupers);  // mirror domainsupers in case view is on device memory
  Kokkos::deep_copy(h_domainsupers, domainsupers);

  Kokkos::parallel_for("initialise_gbxs_on_host", HostTeamPolicy(ngbxs, Kokkos::AUTO()),
                       [=](const HostTeamMember &team_member) {
                         const auto ii = team_member.league_rank();

                         const auto gbx = gen(team_member, ii, gbxmaps, h_domainsupers);

                         /* use 1 thread on host to write gbx to view */
                         team_member.team_barrier();
                         if (team_member.team_rank() == 0) {
                           h_gbxs(ii) = gbx;
                         }
                       });
}

#endif  // LIBS_RUNCLEO_CREATEGBXS_HPP_
