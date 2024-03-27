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
 * Last Modified: Wednesday 27th March 2024
 * Modified By: CB
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
#include <Kokkos_Random.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "./kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "observers2/observers.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

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
 * @tparam Obs Type of the Observer.
 */
template <GridboxMaps GbxMaps, MicrophysicalProcess Microphys, Motion<GbxMaps> M, Observer Obs>
class SDMMethods {
 private:
  unsigned int couplstep; /**< Coupling timestep. */
  MoveSupersInDomain<GbxMaps, M> movesupers;
  /**< object for super-droplets' MoveSupersInDomain with certain type of Motion. */

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
    const auto t_next = (unsigned int)(!(next_mdl < next_move) ? next_move : next_mdl);

    return t_next;
  }

  /* move superdroplets (including movement between
  gridboxes) according to movesupers struct */
  /**
   * @brief Move superdroplets according to the `movesupers` struct.
   *
   * This function moves superdroplets, including their movement between
   * gridboxes, according to the `movesupers` struct. `movesupers` is an
   * instance of the MoveSupersInDomain templated type with a certain instance of a type of
   * GridboxMaps and super-droplets' Motion.
   *
   * @param t_sdm Current timestep for SDM.
   * @param d_gbxs View of gridboxes on device.
   * @param totsupers View of total superdroplets.
   */
  void superdrops_movement(const unsigned int t_sdm, viewd_gbx d_gbxs,
                           const viewd_supers totsupers) const {
    movesupers.run_step(t_sdm, gbxmaps, d_gbxs, totsupers);
  }

 public:
  GbxMaps gbxmaps; /**< object that is type of GridboxMaps. */
  Obs obs;         /**< object that is type of Observer. */

  /**
   * @struct SDMMicrophysics
   * @brief Structure for encapsulating the microphysics process in SDM.
   *
   * The `operator()` is called for SDM microphysics, and it uses Kokkos
   * parallel_for for parallelized execution. struct required so that
   * capture by value KOKKOS_CLASS_LAMBDA (ie. [=] on CPU) only captures
   * objects relevant to microphysics and not other members of SDMMethods
   * (which may not be GPU compatible).
   *
   * @tparam Microphys Type of the MicrophysicalProcess.
   */
  struct SDMMicrophysics {
    Microphys microphys; /**< type of MicrophysicalProcess. */

    /**
     * @brief run SDM microphysics for each gridbox (using sub-timestepping routine).
     *
     * This function runs SDM microphysics for each gridbox using a sub-timestepping routine.
     *  Kokkos::parallel_for is nested parallelism within parallelised loop over gridboxes,
     * serial equivalent is simply: `for (size_t ii(0); ii < ngbxs; ++ii) { [...] }`
     *
     * @param t_sdm Current timestep for SDM.
     * @param t_next Next timestep for SDM.
     * @param d_gbxs View of gridboxes on device.
     * @param genpool Random number generator pool.
     */
    void operator()(const unsigned int t_sdm, const unsigned int t_next, const viewd_gbx d_gbxs,
                    GenRandomPool genpool) const {
      // TODO(all) use scratch space for parallel region?
      const size_t ngbxs(d_gbxs.extent(0));
      Kokkos::parallel_for(
          "sdm_microphysics", TeamPolicy(ngbxs, Kokkos::AUTO()),
          KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
            const int ii = team_member.league_rank();

            auto supers(d_gbxs(ii).supersingbx());
            for (unsigned int subt = t_sdm; subt < t_next; subt = microphys.next_step(subt)) {
              supers = microphys.run_step(team_member, subt, supers, d_gbxs(ii).state, genpool);
            }
          });
    }
  } sdm_microphysics;
  /**< instance of SDMMicrophysics, operator is call of SDM microphysics */

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
             const M movesupers, const Obs obs)
      : couplstep(couplstep),
        movesupers(movesupers),
        gbxmaps(gbxmaps),
        obs(obs),
        sdm_microphysics({microphys}) {}

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
   * @brief Prepare CLEO SDM for timestepping.
   *
   * This function prepares the CLEO SDM for timestepping by
   * calling the `before_timestepping` function of the observer.
   *
   * @param h_gbxs View of gridboxes on host.
   */
  void prepare_to_timestep(const viewh_constgbx h_gbxs) const { obs.before_timestepping(h_gbxs); }

  /**
   * @brief Execute at the start of each coupled model timestep.
   *
   * This function is called at the start of each coupled model timestep
   * (i.e. at start of `t_mdl`) and includes calls to the observer's `at_start_step`
   * function for both the domain and individual gridboxes.
   *
   * @param t_mdl Current timestep of the coupled model.
   * @param d_gbxs View of gridboxes on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs) const {
    const viewd_constsupers totsupers(d_gbxs(0).domain_totsupers_readonly());  // check compatible

    obs.at_start_step(t_mdl, d_gbxs, totsupers);
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
   * @param totsupers View of total superdroplets.
   * @param genpool Random number generator pool.
   */
  void run_step(const unsigned int t_mdl, const unsigned int t_mdl_next, viewd_gbx d_gbxs,
                const viewd_supers totsupers, GenRandomPool genpool) const {
    unsigned int t_sdm(t_mdl);
    while (t_sdm < t_mdl_next) {
      const auto t_sdm_next = next_sdmstep(t_sdm, t_mdl_next);

      superdrops_movement(t_sdm, d_gbxs, totsupers);         // on host and device
      sdm_microphysics(t_sdm, t_sdm_next, d_gbxs, genpool);  // on device

      t_sdm = t_sdm_next;
    }
  }
};

#endif  // LIBS_RUNCLEO_SDMMETHODS_HPP_
