/*
 * ----- CLEO -----
 * File: movesupersindomain.hpp
 * Project: sdmdomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 16th October 2023
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
{
private:
  M motion;

  void move_superdrops_in_domain(const unsigned int t_sdm,
                                 const GridboxMaps &gbxmaps,
                                 viewd_gbx d_gbxs,
                                 viewd_supers supers) const
  {
    // update sd index etc. fro all SDs in gbxs
    motion.update_superdrop_coords(t_sdm);
  }

public:
  MoveSupersInDomain(const M motion)
      : motion(motion) {}

  unsigned int next_step(const unsigned int t_sdm) const
  {
    return motion.next_step(t_sdm);
  }

  void run_step(const unsigned int t_sdm,
                const GridboxMaps &gbxmaps,
                viewd_gbx d_gbxs, 
                viewd_supers supers) const
  {
    if (motion.on_step(t_sdm))
    {
      move_superdrops_in_domain(t_sdm, gbxmaps, d_gbxs, supers);
    }
  };
};

#endif // MOVESUPERSINDOMAIN_HPP