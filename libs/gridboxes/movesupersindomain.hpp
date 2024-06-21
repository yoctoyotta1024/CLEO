/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: movesupersindomain.hpp
 * Project: gridboxes
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality related to moving superdroplets
 * (both updating their spatial coordinates and
 * moving them between gridboxes)
 */

#ifndef LIBS_GRIDBOXES_MOVESUPERSINDOMAIN_HPP_
#define LIBS_GRIDBOXES_MOVESUPERSINDOMAIN_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/sortsupers.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "superdrops/superdrop.hpp"

/*
struct for functionality to move superdroplets throughtout
the domain by updating their spatial coordinates (according to
some type of Motion) and then moving them between gridboxes
after updating their gridbox indexes concordantly
*/
template <GridboxMaps GbxMaps, Motion<GbxMaps> M, typename BoundaryConditions>
struct MoveSupersInDomain {
  /*
  EnactMotion struct encapsulates motion so that parallel loops with KOKKOS_CLASS_LAMBDA
  (ie. [=] on CPUs) functors only captures motion and not other members of MoveSupersInDomain
  coincidentally (which may not be GPU compatible).
  */
  struct EnactMotion {
    M motion;

    /* enact steps (1) and (2) movement of superdroplets for 1 gridbox:
    (1) update their spatial coords according to type of motion. (device)
    (2) update their sdgbxindex accordingly (device).
    Kokkos::parallel_for([...]) is equivalent to:
    for (size_t kk(0); kk < supers.extent(0); ++kk) {[...]}
    when in serial */
    KOKKOS_INLINE_FUNCTION
    void move_supers_in_gbx(const TeamMember &team_member, const unsigned int gbxindex,
                            const GbxMaps &gbxmaps, const State &state,
                            const subviewd_supers supers) const {
      const size_t nsupers(supers.extent(0));
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, nsupers),
                           [&, this](const size_t kk) {
                             /* step (1) */
                             motion.superdrop_coords(gbxindex, gbxmaps, state, supers(kk));

                             /* step (2) */
                             motion.superdrop_gbx(gbxindex, gbxmaps, supers(kk));
                           });
    }

    /* enact steps (1) and (2) movement of superdroplets
    throughout domain (i.e. for all gridboxes):
    (1) update their spatial coords according to type of motion. (device)
    (2) update their sdgbxindex accordingly (device).
    Kokkos::parallel_for([...]) is equivalent to:
    for (size_t ii(0); ii < ngbxs; ++ii) {[...]}
    when in serial */
    void move_supers_in_gridboxes(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs) const {
      const size_t ngbxs(d_gbxs.extent(0));

      Kokkos::parallel_for(
          "move_supers_in_gridboxes", TeamPolicy(ngbxs, Kokkos::AUTO()),
          KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
            const auto ii = team_member.league_rank();

            auto &gbx(d_gbxs(ii));
            move_supers_in_gbx(team_member, gbx.get_gbxindex(), gbxmaps, gbx.state,
                               gbx.supersingbx());
          });
    }

    /* (re)sorting supers based on their gbxindexes and then updating the span for each gridbox
    accordingly.
    Kokkos::parallel_for([...]) (on host) is equivalent to:
    for (size_t ii(0); ii < ngbxs; ++ii){[...]}
    when in serial
    _Note:_ totsupers is view of all superdrops (both in and out of bounds of domain).
    */
    void move_supers_between_gridboxes(const viewd_gbx d_gbxs, const viewd_supers totsupers) const {
      sort_supers(totsupers);

      const size_t ngbxs(d_gbxs.extent(0));
      Kokkos::parallel_for(
          "move_supers_between_gridboxes", TeamPolicy(ngbxs, Kokkos::AUTO()),
          KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
            const auto ii = team_member.league_rank();

            auto &gbx(d_gbxs(ii));
            gbx.supersingbx.set_refs(team_member);
          });

      // /* optional (expensive!) test to raise error if
      // superdrops' gbxindex doesn't match gridbox's gbxindex */
      // for (size_t ii(0); ii < ngbxs; ++ii)
      // {
      //   d_gbxs(ii).supersingbx.iscorrect();
      // }
    }
  } enactmotion;

  MoveSupersInDomain(const M mtn, const BoundaryConditions boundary_conditions)
      : enactmotion({mtn}), apply_domain_boundary_conditions(boundary_conditions) {}

  /* extra constructor useful to help when compiler cannot deduce type of GBxMaps */
  MoveSupersInDomain(const GbxMaps &gbxmaps, const M mtn,
                     const BoundaryConditions boundary_conditions)
      : MoveSupersInDomain(mtn, boundary_conditions) {}

  /* returns time when superdroplet motion is
  next due to occur given current time, t_sdm */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const {
    return enactmotion.motion.next_step(t_sdm);
  }

  /*
   * if current time, t_sdm, is time when superdrop motion should occur, enact movement of
   * superdroplets throughout domain.
   *
   * @param totsupers View of all superdrops (both in and out of bounds of domain).
   *
   */
  void run_step(const unsigned int t_sdm, const GbxMaps &gbxmaps, viewd_gbx d_gbxs,
                const viewd_supers totsupers, const SDMMonitor auto mo) const {
    if (enactmotion.motion.on_step(t_sdm)) {
      move_superdrops_in_domain(t_sdm, gbxmaps, d_gbxs, totsupers);
      mo.monitor_motion(d_gbxs);
    }
  }

 private:
  BoundaryConditions apply_domain_boundary_conditions;

  /* enact movement of superdroplets throughout domain in three stages:
  (1) update their spatial coords according to type of motion. (device)
  (2) update their sdgbxindex accordingly (device)
  (3) move superdroplets between gridboxes (host)
  (4) (optional) apply domain boundary conditions (host and device)
  _Note:_ totsupers is view of all superdrops (both in and out of bounds of domain).
  // TODO(all) use tasking to convert all 3 team policy
  // loops from first two function calls into 1 loop?
  */
  void move_superdrops_in_domain(const unsigned int t_sdm, const GbxMaps &gbxmaps, viewd_gbx d_gbxs,
                                 const viewd_supers totsupers) const {
    /* steps (1 - 2) */
    enactmotion.move_supers_in_gridboxes(gbxmaps, d_gbxs);

    /* step (3) */
    enactmotion.move_supers_between_gridboxes(d_gbxs, totsupers);

    /* step (4) */
    apply_domain_boundary_conditions(gbxmaps, d_gbxs, totsupers);
  }
};

#endif  // LIBS_GRIDBOXES_MOVESUPERSINDOMAIN_HPP_
