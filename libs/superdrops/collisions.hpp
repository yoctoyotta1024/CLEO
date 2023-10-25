/*
 * ----- CLEO -----
 * File: collisions.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for modelling collision
 * microphysical processes in SDM
 * e.g. collision-coalescence
 */

#ifndef COLLISIONS_HPP
#define COLLISIONS_HPP

#include <iostream>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "./microphysicalprocess.hpp"

struct DoCollisions
{
private:
  KOKKOS_FUNCTION void do_collisions(const unsigned int subt) const;

public:

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int subt) const
  /* this operator is used as an "adaptor" for using
  collisions as the MicrophysicsFunction type in a
  ConstTstepMicrophysics instance (*hint* which itself
  satsifies the MicrophysicalProcess concept) */
  {
    do_collisions(subt);
  }

};

inline MicrophysicalProcess auto
Collisions(const unsigned int interval)
/* constructs Microphysical Process for collisions 
of superdroplets with a constant timestep 'interval'
given the "do_collisions" function-like type */
{
  return ConstTstepMicrophysics(interval, DoCollisions{});
}

#endif // COLLISIONS_HPP