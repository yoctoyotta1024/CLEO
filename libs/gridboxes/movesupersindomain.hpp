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
 * Last Modified: Saturday 27th January 2024
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

#include <concepts>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./gridbox.hpp"
#include "./gridboxmaps.hpp"
#include "./sortsupers.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

/* struct for functionality to move superdroplets throughtout
the domain by updating their spatial coordinates (according to
some type of Motion) and then moving them between gridboxes
after updating their gridbox indexes concordantly */
template <GridboxMaps GbxMaps, Motion<GbxMaps> M>
struct MoveSupersInDomain {
  M motion;

  /* enact steps (1) and (2) movement of superdroplets for 1 gridbox:
  (1) update their spatial coords according to type of motion. (device)
  (1b) optional detect precipitation (device)
  (2) update their sdgbxindex accordingly (device).
  Kokkos::parallel_for([...]) is equivalent to:
  for (size_t kk(0); kk < supers.extent(0); ++kk) {[...]}
  when in serial */
  KOKKOS_INLINE_FUNCTION
  void move_supers_in_gbx(const TeamMember &team_member, const unsigned int gbxindex,
                          const GbxMaps &gbxmaps, const State &state,
                          const subviewd_supers supers) const {
    const size_t nsupers(supers.extent(0));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, nsupers), [&, this](const size_t kk) {
      /* step (1) */
      motion.superdrop_coords(gbxindex, gbxmaps, state, supers(kk));

      /* optional step (1b) */
      // gbx.detectors -> detect_precipitation(area, drop); // TODO(CB) detectors

      /* step (2) */
      motion.superdrop_gbx(gbxindex, gbxmaps, supers(kk));
    });
  }

  /* enact steps (1) and (2) movement of superdroplets
  throughout domain (i.e. for all gridboxes):
  (1) update their spatial coords according to type of motion. (device)
  (1b) optional detect precipitation (device)
  (2) update their sdgbxindex accordingly (device).
  Kokkos::parallel_for([...]) is equivalent to:
  for (size_t ii(0); ii < ngbxs; ++ii) {[...]}
  when in serial */
  void move_supers_in_gridboxes(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));

    Kokkos::parallel_for(
        "move_supers_in_gridboxes", TeamPolicy(ngbxs, Kokkos::AUTO()),
        KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
          const int ii = team_member.league_rank();

          auto &gbx(d_gbxs(ii));
          move_supers_in_gbx(team_member, gbx.get_gbxindex(), gbxmaps, gbx.state,
                             gbx.supersingbx());
        });
  }

  /* (re)sorting supers based on their gbxindexes and
  then updating the span for each gridbox accordingly.
  Kokkos::parallel_for([...]) (on host) is equivalent to:
  for (size_t ii(0); ii < ngbxs; ++ii){[...]}
  when in serial */
  void move_supers_between_gridboxes(const viewd_gbx d_gbxs, const viewd_supers totsupers) const {
    sort_supers(totsupers);

    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for(
        "move_supers_between_gridboxes", TeamPolicy(ngbxs, Kokkos::AUTO()),
        KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
          const int ii = team_member.league_rank();

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

  /* enact movement of superdroplets throughout domain in three stages:
  (1) update their spatial coords according to type of motion. (device)
  (1b) optional detect precipitation (device)
  (2) update their sdgbxindex accordingly (device)
  (3) move superdroplets between gridboxes (host) */
  // TODO(all) use tasking to convert all 3 team policy
  // loops in these function calls into 1 loop?
  void move_superdrops_in_domain(const unsigned int t_sdm, const GbxMaps &gbxmaps, viewd_gbx d_gbxs,
                                 const viewd_supers totsupers) const {
    /* steps (1 - 2) */
    move_supers_in_gridboxes(gbxmaps, d_gbxs);

    /* step (3) */
    move_supers_between_gridboxes(d_gbxs, totsupers);
  }

  explicit MoveSupersInDomain(const M i_motion) : motion(i_motion) {}

  /* returns time when superdroplet motion is
  next due to occur given current time, t_sdm */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const { return motion.next_step(t_sdm); }

  /* if current time, t_sdm, is time when superdrop
  motion should occur, enact movement of
  superdroplets throughout domain */
  void run_step(const unsigned int t_sdm, const GbxMaps &gbxmaps, viewd_gbx d_gbxs,
                const viewd_supers totsupers) const {
    if (motion.on_step(t_sdm)) {
      move_superdrops_in_domain(t_sdm, gbxmaps, d_gbxs, totsupers);
    }
  }
};

#endif  // LIBS_GRIDBOXES_MOVESUPERSINDOMAIN_HPP_
