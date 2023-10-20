/*
 * ----- CLEO -----
 * File: sdmmethods.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct wrapping the core ingredients of the Super-droplet Model
 * (SDM) part of CLEO to enact on super-droplets and gridboxes
 */

#ifndef SDMMETHODS_HPP
#define SDMMETHODS_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "../kokkosaliases.hpp"
#include "./coupleddynamics.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "observers/observers.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

template <CoupledDynamics CD, GridboxMaps GbxMaps,
          MicrophysicalProcess Microphys,
          Motion M, Observer Obs>
struct SDMMethods
{
private:
  unsigned int couplstep; // coupled timestep

  KOKKOS_INLINE_FUNCTION
  unsigned int next_sdmstep(const unsigned int t_sdm,
                            const unsigned int next_mdl) const
  /* given current timestep, t_sdm, work out which event
  (motion or one complete step) is next to occur and return
  the time of the sooner event, (ie. next t_move or t_mdl) */
  {
    const unsigned int next_move(movesupers.next_step(t_sdm));        // t of next sdm movement
    const unsigned int t_next(!(next_mdl <
                                next_move)
                                  ? next_move
                                  : next_mdl);

    return t_next; // return smaller of two unsigned ints (see std::min)
  }

  KOKKOS_INLINE_FUNCTION
  void superdrops_movement(const unsigned int t_sdm,
                           const viewd_gbx d_gbxs,
                           const viewd_supers supers) const
  /* move superdroplets (including movement between
  gridboxes) according to movesupers struct */
  {
    movesupers.run_step(t_sdm, gbxmaps, d_gbxs, supers);
  }

  KOKKOS_INLINE_FUNCTION
  void sdm_microphysics(const unsigned int t_sdm,
                        const unsigned int t_next,
                        const viewd_gbx d_gbxs) const
  /* enact SDM microphysics for each gridbox
  (using sub-timestepping routine) */
  {
    // loop over gbxs for(gbx : gbxs)
    {
      for (unsigned int subt = t_sdm; subt < t_next;
           subt = microphys.next_step(subt))
      {
        microphys.run_step(subt);
      }
    }
  }

public:
  GbxMaps gbxmaps;                  // maps from gridbox indexes to domain coordinates
  Microphys microphys;              // microphysical process
  MoveSupersInDomain<M> movesupers; // super-droplets' motion in domain
  Obs obs;                          // observer

  SDMMethods(const CD &coupldyn,
             const GbxMaps gbxmaps,
             const Microphys microphys,
             const MoveSupersInDomain<M> movesupers,
             const Obs obs)
      : couplstep(coupldyn.get_couplstep()),
        gbxmaps(gbxmaps),
        microphys(microphys),
        movesupers(movesupers),
        obs(obs) {}

  KOKKOS_INLINE_FUNCTION
  auto get_couplstep() const { return couplstep; }

  void prepare_to_timestep(const viewh_constgbx h_gbxs) const
  /* prepare CLEO SDM for timestepping */
  {
    obs.before_timestepping(h_gbxs);
  }

  void receive_dynamics(const CD &coupldyn,
                        const viewh_gbx h_gbxs) const {}
  /* update Gridboxes' states using information
  received from coupldyn */

  void send_dynamics(const CD &coupldyn,
                     const viewh_constgbx h_gbxs) const {}
  /* send information from Gridboxes' states to coupldyn */

  KOKKOS_INLINE_FUNCTION
  void run_step(const unsigned int t_mdl,
                const unsigned int t_mdl_next,
                const viewd_gbx d_gbxs,
                const viewd_supers supers) const
  /* run CLEO SDM (on device) from time t_mdl to t_mdl_next
  with sub-timestepping routine for super-droplets'
  movement and microphysics */
  {
    unsigned int t_sdm(t_mdl);
    while (t_sdm < t_mdl_next)
    {
      unsigned int t_sdm_next(next_sdmstep(t_sdm, t_mdl_next));

      superdrops_movement(t_sdm, d_gbxs, supers);
      sdm_microphysics(t_sdm, t_sdm_next, d_gbxs);

      t_sdm = t_sdm_next;
    }
  }
};

#endif // SDMMETHODS_HPP