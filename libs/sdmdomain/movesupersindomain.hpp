/*
 * ----- CLEO -----
 * File: movesupersindomain.hpp
 * Project: sdmdomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
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

#include <map>
#include <utility>
#include <stdexcept>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "./gridbox.hpp"
#include "./gridboxmaps.hpp"
#include "superdrops/superdrop.hpp"

struct SuperdropMotion
{
private:
  const unsigned int interval;

public:
  SuperdropMotion(const unsigned int motionstep)
      : interval(motionstep) {}

  bool on_step(const unsigned int t_sdm) const
  {
    return t_sdm % interval == 0;
  }

  unsigned int next_step(const unsigned int t_sdm) const
  {
    return t_sdm + interval;
  }
};

struct MoveSupersInDomain
{
private:
  SuperdropMotion sdmotion;

  void move_supers_in_domain(const unsigned int t_sdm,
                             const GridboxMaps &gbxmaps,
                             Gridboxes &gbxs,
                             Superdrops &supers) const
  {
    std::cout << "move @ t = " << t_sdm << "\n";
  }

public:
  MoveSupersInDomain(const unsigned int motionstep)
      : sdmotion(motionstep) {}

  unsigned int next_step(const unsigned int t_sdm) const
  {
    return sdmotion.next_step(t_sdm);
  }

  void run_step(const unsigned int t_sdm,
                const GridboxMaps &gbxmaps,
                Gridboxes &gbxs,
                Superdrops &supers) const
  {
    if (sdmotion.on_step(t_sdm))
    {
      move_supers_in_domain(t_sdm, gbxmaps, gbxs, supers);
    }
  };
};

#endif // MOVESUPERSINDOMAIN_HPP