/*
 * ----- CLEO -----
 * File: movesupersindomain.hpp
 * Project: gridboxes
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 7th November 2023
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
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

template <Motion M>
struct MoveSupersInDomain
/* struct for functionality to move superdroplets throughtout
the domain by updating their spatial coordinates (according to 
some type of Motion) and then moving them between gridboxes 
after updating their gridbox indexes concordantly */
{
private:
  M motion;

  KOKKOS_INLINE_FUNCTION
  unsigned int update_superdrop_gbxindex() const
  {
    // TODO (put into seperate templated struct too)
    return 0;
  }

  void move_superdroplets_between_gridboxes() const
  {
    // TODO (put into seperate templated struct too)
  }

  void move_superdrops_in_domain(const unsigned int t_sdm,
                                 const GridboxMaps auto &gbxmaps,
                                 const viewd_gbx d_gbxs,
                                 const viewd_supers totsupers) const
  /* enact movement of superdroplets throughout domain in three stages:
  (1) update their spatial coords according to type of motion.
  (1b) optional detect precipitation
  (2) update their sdgbxindex accordingly
  (3) move superdroplets between gridboxes */
  {
    const size_t ngbxs(d_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      const subviewd_supers supers(d_gbxs(ii).supersingbx());
      for (size_t kk(0); kk < supers.extent(0); ++kk)
      {
        Superdrop &drop(supers(kk));

        /* step (1) */
        motion.update_superdrop_coords();

        /* optional step (1b) */
        // gbx.detectors -> detect_precipitation(area, drop); // TODO (detectors)

        /* step (2) */
        drop.set_sdgbxindex(update_superdrop_gbxindex());
      }
    }

    /* step (3) */
    move_superdroplets_between_gridboxes();
  }

public:
  MoveSupersInDomain(const M motion)
      : motion(motion) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const
  /* returns time when superdroplet motion is
  next due to occur given current time, t_sdm */
  {
    return motion.next_step(t_sdm);
  }

  void run_step(const unsigned int t_sdm,
                const GridboxMaps auto &gbxmaps,
                const viewd_gbx d_gbxs,
                const viewd_supers totsupers) const
  /* if current time, t_sdm, is time when superdrop
  motion should occur, enact movement of
  superdroplets throughout domain */
  {
    if (motion.on_step(t_sdm))
    {
      move_superdrops_in_domain(t_sdm, gbxmaps, d_gbxs, totsupers);
    }
  };
};

#endif // MOVESUPERSINDOMAIN_HPP