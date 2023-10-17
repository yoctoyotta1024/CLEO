/*
 * ----- CLEO -----
 * File: movesupersindomain.hpp
 * Project: gridboxes
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
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

  void move_superdrops_in_domain(const unsigned int t_sdm,
                                 const GridboxMaps auto &gbxmaps,
                                 viewd_gbx d_gbxs,
                                 viewd_supers supers) const
  /* enact movement of superdroplets throughout domain in three stages:
  1) update their spatial coords according to type of motion.
  2) update their sd_gbxindex accordingly
  3) move superdroplets between gridboxes */
  {
    // update sd index etc. fro all SDs in gbxs
    motion.update_superdrop_coords(t_sdm);
  }

public:
  MoveSupersInDomain(const M motion)
      : motion(motion) {}

  unsigned int next_step(const unsigned int t_sdm) const
  /* returns time when superdroplet motion is
  next due to occur given current time, t_sdm */
  {
    return motion.next_step(t_sdm);
  }

  void run_step(const unsigned int t_sdm,
                const GridboxMaps auto &gbxmaps,
                viewd_gbx d_gbxs,
                viewd_supers supers) const
  /* if current time, t_sdm, is time when superdrop
  motion should occur, enact movement of
  superdroplets throughout domain */
  {
    if (motion.on_step(t_sdm))
    {
      move_superdrops_in_domain(t_sdm, gbxmaps, d_gbxs, supers);
    }
  };
};

#endif // MOVESUPERSINDOMAIN_HPP