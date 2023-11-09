/*
 * ----- CLEO -----
 * File: coalescence.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 9th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * class that enacts collision-coalescence events
 * in superdroplet model according to Shima et al. 2009.
 * Coalescence struct satisfies SDinGBxPairEnactX concept
 * used in Collisions struct */

#ifndef COALESCENCE_HPP
#define COALESCENCE_HPP

#include <Kokkos_Core.hpp>

#include "./collisions.hpp"
#include "./kokkosaliases_sd.hpp"
#include "./microphysicalprocess.hpp"
#include "./superdrop.hpp"
struct DoCoalescence
{
private:
public:
};

inline MicrophysicalProcess auto
Collisions(const unsigned int interval, const Kernel kernel)
/* constructs Microphysical Process for collisions 
of superdroplets with a constant timestep 'interval'
given the "do_collisions" function-like type */
{
  const DoCollisions<Kernel, DoCoalescence> colls(kernel,
                                                  DoCoalescence{}); // TODO use actualy kernel 
  return ConstTstepMicrophysics(interval, colls);
}

#endif // COLLISIONCOALESCENCE_HPP