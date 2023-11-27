/*
 * ----- CLEO -----
 * File: movesupersindomain.hpp
 * Project: gridboxes
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 8th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Functionality related to moving superdroplets
 * (both updating their spatial coordinates and
 * moving them between gridboxes)
 */

#ifndef MOVESUPERSINDOMAIN_HPP
#define MOVESUPERSINDOMAIN_HPP

#include <concepts>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./gridbox.hpp"
#include "./gridboxmaps.hpp"
#include "./sortsupers.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

template <GridboxMaps GbxMaps,
          Motion<GbxMaps> M>
struct MoveSupersInDomain
/* struct for functionality to move superdroplets throughtout
the domain by updating their spatial coordinates (according to
some type of Motion) and then moving them between gridboxes
after updating their gridbox indexes concordantly */
{
  M motion;

  KOKKOS_INLINE_FUNCTION
  void move_supers_in_gbx(const TeamMember &team_member,
                          const unsigned int gbxindex,
                          const GbxMaps &gbxmaps,
                          const State &state,
                          const subviewd_supers supers) const
  /* enact steps (1) and (2) movement of superdroplets for 1 gridbox:
  (1) update their spatial coords according to type of motion. (device)
  (1b) optional detect precipitation (device)
  (2) update their sdgbxindex accordingly (device).
  Kokkos::parallel_for is equivalent to:
  for (size_t kk(0); kk < supers.extent(0); ++kk) {[...]} in serial */
  {
    const size_t nsupers(supers.extent(0));
    Kokkos::parallel_for(
        "move_superdrops_in_gbx",
        Kokkos::TeamThreadRange(team_member, nsupers),
        [=](const size_t kk)
        {
          /* step (1) */
          motion.update_superdrop_coords(gbxindex,
                                         gbxmaps,
                                         state,
                                         supers(kk));

          /* optional step (1b) */
          // gbx.detectors -> detect_precipitation(area, drop); // TODO detectors

          /* step (2) */
          motion.update_superdrop_gbxindex(gbxindex,
                                           gbxmaps,
                                           supers(kk));
        });
  }

  void move_supers_in_gridboxes(const GbxMaps &gbxmaps,
                                const viewd_gbx d_gbxs) const
  /* enact steps (1) and (2) movement of superdroplets
  throughout domain (i.e. for all gridboxes):
  (1) update their spatial coords according to type of motion. (device)
  (1b) optional detect precipitation (device)
  (2) update their sdgbxindex accordingly (device).
  Kokkos::parallel_for is equivalent to:
  for (size_t ii(0); ii < ngbxs; ++ii) {[...]} in serial */
  {
    const size_t ngbxs(d_gbxs.extent(0));

    Kokkos::parallel_for(
        "move_superdrops_in_domain",
        TeamPolicy(ngbxs, Kokkos::AUTO()),
        KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
          const int ii = team_member.league_rank();

          auto &gbx(d_gbxs(ii));
          move_supers_in_gbx(team_member,
                             gbx.get_gbxindex(),
                             gbxmaps,
                             gbx.state,
                             gbx.supersingbx());
        });
  }

  void move_supers_between_gridboxes(const viewh_gbx h_gbxs,
                                     const viewd_supers totsupers) const
  /* (re)sorting supers based on their gbxindexes and
  then updating the span for each gridbox accordingly */
  {
    sort_supers(totsupers);

    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii) // TODO parallelise on host?
    {
      h_gbxs(ii).supersingbx.set_refs();
      // h_gbxs(ii).supersingbx.iscorrect(); // (expensive!) optional test to raise error if superdrops' gbxindex doesn't match gridbox's gbxindex
    }
  }

  void move_superdrops_in_domain(const unsigned int t_sdm,
                                 const GbxMaps &gbxmaps,
                                 dualview_gbx gbxs,
                                 const viewd_supers totsupers) const
  /* enact movement of superdroplets throughout domain in three stages:
  (1) update their spatial coords according to type of motion. (device)
  (1b) optional detect precipitation (device)
  (2) update their sdgbxindex accordingly (device)
  (3) move superdroplets between gridboxes (host) */
  {
    /* steps (1 - 2) */
    gbxs.sync_device(); // get device up to date with host
    move_supers_in_gridboxes(gbxmaps, gbxs.view_device());
    gbxs.modify_device(); // mark device view of gbxs as modified

    /* step (3) */
    gbxs.sync_host(); // get device up to date with host
    move_supers_between_gridboxes(gbxs.view_host(), totsupers);
    gbxs.modify_host(); // mark device view of gbxs as modified

    gbxs.sync_device(); // get device up to date with host
  }

  MoveSupersInDomain(const M i_motion)
      : motion(i_motion) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const
  /* returns time when superdroplet motion is
  next due to occur given current time, t_sdm */
  {
    return motion.next_step(t_sdm);
  }

  void run_step(const unsigned int t_sdm,
                const GbxMaps &gbxmaps,
                dualview_gbx gbxs,
                const viewd_supers totsupers) const
  /* if current time, t_sdm, is time when superdrop
  motion should occur, enact movement of
  superdroplets throughout domain */
  {
    if (motion.on_step(t_sdm))
    {
      move_superdrops_in_domain(t_sdm, gbxmaps, gbxs, totsupers);
    }
  };
};

#endif // MOVESUPERSINDOMAIN_HPP