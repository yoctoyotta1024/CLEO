/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: sdmmethods.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 27th January 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct wrapping the core ingredients of the Super-droplet Model
 * (SDM) part of CLEO to enact on super-droplets and gridboxes
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
#include "observers/observers.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

/**
 * @class SDMMethods
 * @brief Struct wrapping the core ingredients of the Super-droplet Model (SDM) part of CLEO.
 *
 * This struct encapsulates the essential components of the Super-droplet Model (SDM) in the CLEO coupled model.
 * It includes components for handling gridboxes, super-droplets' motion, microphysics, and observers.
 *
 * @tparam GbxMaps Type of the GridboxMaps.
 * @tparam Microphys Type of the MicrophysicalProcess.
 * @tparam M Type of the Motion.
 * @tparam Obs Type of the Observer.
 */
template <GridboxMaps GbxMaps, MicrophysicalProcess Microphys, Motion<GbxMaps> M, Observer Obs>
class SDMMethods {
 private:
  unsigned int couplstep;                     // coupled timestep
  MoveSupersInDomain<GbxMaps, M> movesupers;  // super-droplets' motion in domain

  /* given current timestep, t_sdm, work out which event
  (motion or one complete step) is next to occur and return
  the time of the sooner event, (ie. next t_move or t_mdl) */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_sdmstep(const unsigned int t_sdm, const unsigned int next_mdl) const {
    const auto next_move = (unsigned int)movesupers.next_step(t_sdm);  // t of next sdm movement
    const auto t_next = (unsigned int)(!(next_mdl < next_move) ? next_move : next_mdl);

    return t_next;  // return smaller of two unsigned ints (see std::min)
  }

  /* move superdroplets (including movement between
  gridboxes) according to movesupers struct */
  void superdrops_movement(const unsigned int t_sdm, viewd_gbx d_gbxs,
                           const viewd_supers totsupers) const {
    movesupers.run_step(t_sdm, gbxmaps, d_gbxs, totsupers);
  }

 public:
  GbxMaps gbxmaps;  // maps from gridbox indexes to domain coordinates
  Obs obs;          // observer

  /* operator is call to microphysics 'sdm_microphysics'. struct
  created so that capture by value KOKKOS_CLASS_LAMBDA
  (ie. [=] on CPU) only captures microphysics and not other
  members of SDMMethods.
  ()*/
  struct SDMMicrophysics {
    Microphys microphys;  // microphysical process

    /* enact SDM microphysics for each gridbox
    (using sub-timestepping routine). Kokkos::parallel for is
    nested parallelism within parallelised loop over gridboxes,
    serial equivalent is simply:
    for (size_t ii(0); ii < ngbxs; ++ii) { [...] } */
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
  } sdm_microphysics;  // operator is call for SDM microphysics

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
   *
   */
  KOKKOS_INLINE_FUNCTION
  auto get_couplstep() const { return couplstep; }

  /* prepare CLEO SDM for timestepping */
  void prepare_to_timestep(const viewh_constgbx h_gbxs) const { obs.before_timestepping(h_gbxs); }

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs) const {
    const viewd_constsupers totsupers(h_gbxs(0).domain_totsupers_readonly());

    obs.at_start_step(t_mdl, h_gbxs, totsupers);

    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii) {
      obs.at_start_step(t_mdl, h_gbxs(ii));
    }
  }

  /* run CLEO SDM (on device) from time t_mdl to t_mdl_next
  with sub-timestepping routine for super-droplets'
  movement and microphysics */
  void run_step(const unsigned int t_mdl, const unsigned int t_mdl_next, viewd_gbx d_gbxs,
                const viewd_supers totsupers, GenRandomPool genpool) const {
    unsigned int t_sdm(t_mdl);
    while (t_sdm < t_mdl_next) {
      unsigned int t_sdm_next(next_sdmstep(t_sdm, t_mdl_next));

      superdrops_movement(t_sdm, d_gbxs, totsupers);         // on host and device
      sdm_microphysics(t_sdm, t_sdm_next, d_gbxs, genpool);  // on device

      t_sdm = t_sdm_next;
    }
  }
};

#endif  // LIBS_RUNCLEO_SDMMETHODS_HPP_
