/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: sdmmethods.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors: Tobias Kölling (TK)
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct wrapping the core ingredients of CLEO's Super-droplet Model (SDM)
 * the microphysical process, motion etc. to enact on super-droplets and gridboxes
 */

#ifndef LIBS_RUNCLEO_SDMMETHODS_HPP_
#define LIBS_RUNCLEO_SDMMETHODS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "./kokkosaliases.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "gridboxes/supersindomain.hpp"
#include "gridboxes/transport_across_domain.hpp"
#include "observers/observers.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "superdrops/superdrop.hpp"

namespace KCS = KokkosCleoSettings;

/**
 * @struct SDMMicrophysicsFunctor
 * @brief Structure for encapsulating the microphysics process in SDM.
 *
 * The `operator()` is called for SDM microphysics, and it uses Kokkos parallel_reduce(...)
 * for parallelized execution over gridboxes and/or superdroplets. Struct ensures parallel region
 * only captures objects relevant to microphysics and not other members of SDMMethods
 * (which may not be GPU compatible).
 *
 * @tparam Microphys Type of the MicrophysicalProcess.
 */
template <MicrophysicalProcess Microphys, SDMMonitor SDMMo>
struct SDMMicrophysicsFunctor {
  const Microphys microphys;          /**< object that is type of MicrophysicalProcess. */
  const unsigned int t_sdm;           /**< current timestep for SDM. */
  const unsigned int t_next;          /**< next timestep for SDM. */
  const viewd_gbx d_gbxs;             /** view of gridboxes on device. */
  const subviewd_supers domainsupers; /**view on device of all superdroplets in all gridboxes. */
  const SDMMo mo;                     /**< object that is type of SDMMonitor to use. */

  /** Note(!): number of superdroplets in supers may change in call to run microphysics and, if so,
   * then it is REQUIRED to sort allsupers outside this functor (i.e. after gbxs parallel loop)
   * and then call set_refs for all gbxs. See sdm_microphysics(...) below. This can occur e.g.
   * due to null superdroplet after collision coalescence of two xi=1 superdroplets.
   *
   * @return true if number of superdroplets changes, false otherwise
   */
  KOKKOS_INLINE_FUNCTION bool run_subtimestepping(const TeamMember& team_member, State& state,
                                                  subviewd_supers supers) const {
    const auto nsupers_before = static_cast<size_t>(supers.extent(0));
    for (unsigned int subt = t_sdm; subt < t_next; subt = microphys.next_step(subt)) {
      supers = microphys.run_step(team_member, subt, supers, state, mo);
    }
    mo.monitor_microphysics(team_member, supers);
    const auto nsupers_after = static_cast<size_t>(supers.extent(0));

    if (nsupers_before == nsupers_after) {
      return false;
    } else {
      return true;  // number of superdroplets has changed
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team_member,
                                         bool& any_nsupers_change) const {
    const auto ii = team_member.league_rank();
    const auto nsupers_change =
        run_subtimestepping(team_member, d_gbxs(ii).state, d_gbxs(ii).supersingbx(domainsupers));
    any_nsupers_change = any_nsupers_change || nsupers_change;
  }
};

/**
 * @class SDMMethods
 * @brief Struct wrapping the core ingredients of the Super-droplet Model (SDM) part of CLEO.
 *
 * This struct encapsulates the essential components of the Super-droplet Model (SDM)
 * in the CLEO coupled model. It includes components for handling gridboxes,
 * super-droplets' motion, microphysics, and observers.
 *
 * @tparam GbxMaps Type of the GridboxMaps.
 * @tparam Microphys Type of the MicrophysicalProcess.
 * @tparam M Type of super-droplets' Motion.
 * @tparam TransportAcrossDomain Type of super-droplets transport across domain.
 * @tparam BoundaryConditions Type of boundary conditions for superdroplet motion
 * @tparam Obs Type of the Observer.
 */
template <GridboxMaps GbxMaps, MicrophysicalProcess Microphys, Motion<GbxMaps> M,
          TransportAcrossDomain<GbxMaps> T, BoundaryConditions<GbxMaps> BCs, Observer Obs>
class SDMMethods {
 private:
  unsigned int couplstep; /**< Coupling timestep. */
  MoveSupersInDomain<GbxMaps, M, T, BCs> movesupers;
  /**< object for super-droplets' MoveSupersInDomain with certain type of Motion, transport and
   * boundary conditions. */

  /**
   * @brief Get the next timestep for SDM.
   *
   * Given the current timestep for SDM (`t_sdm`) and the next timestep for the
   * coupled model (`next_mdl`), this function determines which event
   * (motion or one complete step) will be the next to occur and returns the
   * time of the sooner event (i.e., next `t_move` or `t_mdl`).
   *
   * @param t_sdm Current timestep of SDM.
   * @param next_mdl Next timestep of the coupled model.
   * @return The timestep of the sooner event.
   */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_sdmstep(const unsigned int t_sdm, const unsigned int next_mdl) const {
    const auto next_move = (unsigned int)movesupers.next_step(t_sdm);

    /* return smaller of two unsigned ints (see std::min) */
    const auto t_next = Kokkos::min(next_mdl, next_move);

    return t_next;
  }

  /**
   * @brief Move superdroplets according to the `movesupers` struct.
   *
   * This function moves superdroplets, including their movement between
   * gridboxes and boundary conditions, according to the `movesupers` struct.
   * `movesupers` is an instance of the MoveSupersInDomain templated type with a certain
   * instance of a type of GridboxMaps, super-droplets' Motion, transport and boundary conditions.
   *
   * Kokkos::Profiling are null pointers unless a Kokkos profiler library has been
   * exported to "KOKKOS_TOOLS_LIBS" prior to runtime so the lib gets dynamically loaded.
   *
   * @param t_sdm Current timestep for SDM.
   * @param d_gbxs View of gridboxes on device.
   * @param allsupers View of all superdrops (both in and out of bounds of domain).
   * @param mo Monitor of SDM processes.
   */
  void superdrops_movement(const unsigned int t_sdm, viewd_gbx d_gbxs, SupersInDomain& allsupers,
                           const SDMMonitor auto mo) const {
    Kokkos::Profiling::ScopedRegion region("timestep_sdm_movement");

    allsupers = movesupers.run_step(t_sdm, gbxmaps, d_gbxs, allsupers, mo);
  }

 public:
  GbxMaps gbxmaps;     /**< object that is type of GridboxMaps. */
  Obs obs;             /**< object that is type of Observer. */
  Microphys microphys; /**< object that is type of MicrophysicalProcess. */

  /**
   * @brief run SDM microphysics for each gridbox (using sub-timestepping routine).
   *
   * This function runs SDM microphysics for each gridbox using a sub-timestepping routine.
   * Kokkos::parallel_reduce is nested parallelism within parallelised loop over gridboxes,
   * serial equivalent is simply: `for (size_t ii(0); ii < ngbxs; ++ii) { [...] }` which then
   * returns combined result of "logical or" operation over all gridboxes.
   *
   * @param t_sdm Current timestep for SDM.
   * @param t_next Next timestep for SDM.
   * @param d_gbxs View of gridboxes on device.
   * @param domainsupers View on device of all the superdroplets related to the gridboxes.
   * @param mo SDMMonitor to use.
   * @return true if number of superdroplets in any gridbox has changed during microphysics
   */
  template <SDMMonitor SDMMo>
  bool sdm_microphysics(const unsigned int t_sdm, const unsigned int t_next, const viewd_gbx d_gbxs,
                        const subviewd_supers domainsupers, const SDMMo mo) const {
    // TODO(ALL) use scratch space for parallel region(?)
    const size_t ngbxs(d_gbxs.extent(0));
    const auto functor = SDMMicrophysicsFunctor{microphys, t_sdm, t_next, d_gbxs, domainsupers, mo};

    auto any_nsupers_change = bool{false};
    Kokkos::parallel_reduce("sdm_microphysics", TeamPolicy(ngbxs, KCS::team_size), functor,
                            Kokkos::LOr<bool>(any_nsupers_change));
    return any_nsupers_change;
  }

  /**
   * @brief run SDM microphysics for each gridbox (using sub-timestepping routine).
   *
   * This function is a wrapper around the function which runs SDM microphysics and then will
   * sort all the superdroplets and set the refs of all the gridboxes if the change in number of
   * superdroplets boolean returned from running microphysics is true
   *
   * Kokkos::Profiling are null pointers unless a Kokkos profiler library has been
   * exported to "KOKKOS_TOOLS_LIBS" prior to runtime so the lib gets dynamically loaded.
   *
   * @param t_sdm Current timestep for SDM.
   * @param t_next Next timestep for SDM.
   * @param d_gbxs View of gridboxes on device.
   * @param allsupers Struct to handle superdroplets in the domain.
   * @param mo SDMMonitor to use.
   */
  template <SDMMonitor SDMMo>
  void sdm_microphysics(const unsigned int t_sdm, const unsigned int t_next, const viewd_gbx d_gbxs,
                        SupersInDomain& allsupers, const SDMMo mo) const {
    Kokkos::Profiling::ScopedRegion region("timestep_sdm_microphysics");

    const auto domainsupers = allsupers.domain_supers();
    const auto any_nsupers_change = sdm_microphysics(t_sdm, t_next, d_gbxs, domainsupers, mo);

    if (any_nsupers_change) {
      allsupers.sort_totsupers(d_gbxs);
      const size_t ngbxs(d_gbxs.extent(0));
      Kokkos::parallel_for(
          "microphysics_set_gridboxes_refs", Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
          KOKKOS_LAMBDA(const size_t ii) {
            d_gbxs(ii).supersingbx.set_refs(allsupers.domain_supers());
          });
    }
  }

  /**
   * @brief Constructor for SDMMethods.
   *
   * Initializes SDMMethods with the provided coupling timestep, gridbox maps,
   * microphysics, motion, and observer.
   *
   * @param couplstep Coupling timestep.
   * @param gbxmaps object that is type of GridboxMaps.
   * @param microphys object that is type of MicrophysicalProcess.
   * @param movesupers object that is type of super-droplets' Motion.
   * @param obs object that is type of Observer.
   */
  SDMMethods(const unsigned int couplstep, const GbxMaps gbxmaps, const Microphys microphys,
             const MoveSupersInDomain<GbxMaps, M, T, BCs> movesupers, const Obs obs)
      : couplstep(couplstep),
        movesupers(movesupers),
        gbxmaps(gbxmaps),
        obs(obs),
        microphys(microphys) {}

  /**
   * @brief Get the coupling step value.
   *
   * This function retrieves and returns the size of the coupling timestep.
   *
   * @return The coupling timestep value.
   */
  KOKKOS_INLINE_FUNCTION
  auto get_couplstep() const { return couplstep; }

  /**
   * @brief Given current timestep return time of next coupling event.
   *
   * This function returns the time of the next coupling event.
   *
   * @return The next coupling timestep.
   */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_couplstep(const unsigned int t_mdl) const {
    return ((t_mdl / couplstep) + 1) * couplstep;
  }

  /**
   * @brief Prepare CLEO SDM for timestepping.
   *
   * This function prepares the CLEO SDM for timestepping by
   * calling the `before_timestepping` function of the observer.
   *
   * @param gbxs DualView of gridboxes.
   * @param allsupers Superdroplets in and out of the domain.
   */
  void prepare_to_timestep(const dualview_gbx gbxs, const SupersInDomain& allsupers) const {
    const auto d_gbxs = gbxs.view_device();
    const auto domainsupers = allsupers.domain_supers_readonly();
    obs.before_timestepping(d_gbxs, domainsupers);
  }

  /**
   * @brief Execute at the start of each coupled model timestep.
   *
   * This function is called at the start of each coupled model timestep
   * (i.e. at start of `t_mdl`) and includes calls to the observer's `at_start_step`
   * function for both the domain and individual gridboxes.
   *
   * @param t_mdl Current timestep of the coupled model.
   * @param gbxs Dualview of gridboxes (on host and on device).
   * @param allsupers View of all (inside and outside of domain) superdroplets on device.
   */
  void at_start_step(const unsigned int t_mdl, const dualview_gbx gbxs,
                     const SupersInDomain& allsupers) const {
    const auto d_gbxs = gbxs.view_device();
    const auto domainsupers = allsupers.domain_supers_readonly();
    obs.at_start_step(t_mdl, d_gbxs, domainsupers);
  }

  /**
   * @brief Run CLEO SDM for a specified timestep range.
   *
   * This function runs CLEO SDM on the device from time `t_mdl` to `t_mdl_next`,
   * with a sub-timestepping routine for the super-droplets' movement
   * and microphysics.
   *
   * @param t_mdl Current timestep of the coupled model.
   * @param t_mdl_next Next timestep of the coupled model.
   * @param d_gbxs View of gridboxes on device.
   * @param allsupers View of all superdrops (both in and out of bounds of domain).
   */
  void run_step(const unsigned int t_mdl, const unsigned int t_mdl_next, viewd_gbx d_gbxs,
                SupersInDomain& allsupers) const {
    const SDMMonitor auto mo = obs.get_sdmmonitor();

    unsigned int t_sdm(t_mdl);
    while (t_sdm < t_mdl_next) {
      const auto t_sdm_next = next_sdmstep(t_sdm, t_mdl_next);

      superdrops_movement(t_sdm, d_gbxs, allsupers, mo);           // on host and device
      sdm_microphysics(t_sdm, t_sdm_next, d_gbxs, allsupers, mo);  // on device

      t_sdm = t_sdm_next;
    }
  }
};

#endif  // LIBS_RUNCLEO_SDMMETHODS_HPP_
