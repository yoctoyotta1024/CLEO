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

template <typename S, typename GbxMaps>
concept UpdateSdgbxindex = requires(S s, const unsigned int u,
                                    const GbxMaps &gbxmaps,
                                    Superdrop &drop)
/* concept for all (function-like) types (ie. types
that can be called with some arguments) that can be
called by MoveSupersInDomain for the
"update_superdrop_gbxindex" function (see below) */
{
  {
    s(u, gbxmaps, drop)
  } -> std::same_as<void>;
};

template <GridboxMaps GbxMaps,
          Motion<GbxMaps> M,
          UpdateSdgbxindex<GbxMaps> S>
struct MoveSupersInDomain
/* struct for functionality to move superdroplets throughtout
the domain by updating their spatial coordinates (according to
some type of Motion) and then moving them between gridboxes
after updating their gridbox indexes concordantly */
{
  M motion;
  S update_superdrop_gbxindex;

  void move_supers_between_gridboxes(const viewh_gbx h_gbxs,
                                     const viewd_supers totsupers) const
  /* (re)sorting supers based on their gbxindexes and
  then updating the span for each gridbox accordingly */
  {
    sort_supers(totsupers);

    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      h_gbxs(ii).supersingbx.set_refs();
      // h_gbxs(ii).supersingbx.iscorrect(); // (expensive!) optional test to raise error if superdrops' gbxindex doesn't match gridbox's gbxindex
    }
  }

  void move_supers_in_gridboxes(const GbxMaps &gbxmaps,
                                const viewd_gbx d_gbxs) const
  /* enact movement of superdroplets throughout domain in three stages:
  (1) update their spatial coords according to type of motion. (device)
  (1b) optional detect precipitation (device)
  (2) update their sdgbxindex accordingly (device) */
  {
    const size_t ngbxs(d_gbxs.extent(0));

    /* parallelised version of:
    for (size_t ii(0); ii < ngbxs; ++ii) {[...]} */
    Kokkos::parallel_for(
        "move_superdrops_in_domain",
        Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
        KOKKOS_CLASS_LAMBDA(const size_t ii) {
          const subviewd_supers supers(d_gbxs(ii).supersingbx());
          for (size_t kk(0); kk < supers.extent(0); ++kk)
          {
            const unsigned int gbxindex(d_gbxs(ii).get_gbxindex());
            
            /* step (1) */
            motion.update_superdrop_coords(gbxindex, gbxmaps,
                                           d_gbxs(ii).state,
                                           supers(kk));

            /* optional step (1b) */
            // gbx.detectors -> detect_precipitation(area, drop); // TODO (detectors)

            /* step (2) */
            update_superdrop_gbxindex(gbxindex, gbxmaps, supers(kk));
            // supers(kk).set_sdgbxindex(update_superdrop_gbxindex()); // TODO fill in update func
          }
        });
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

  MoveSupersInDomain(const M i_motion, const S i_sdgbxfunc)
      : motion(i_motion),
        update_superdrop_gbxindex(i_sdgbxfunc) {}

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