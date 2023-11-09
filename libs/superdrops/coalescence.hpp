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
  KOKKOS_FUNCTION
  void operator()(SuperdropWithGbxindex &SDinGBx1,
                  SuperdropWithGbxindex &SDinGBx2,
                  const double prob, const double phi) const
  /* this operator is used as an "adaptor" for using
  DoCoalescence as a function in DoCollisions that
  satistfies the PairEnactX concept */
  {
    // TODO
  }
};

template <PairProbability Probability>
inline MicrophysicalProcess auto
CollCoal(const unsigned int interval,
         const std::function<double(int)> int2realtime,
         const Probability collcoalprob)
/* constructs Microphysical Process for collision-coalescence
of superdroplets with a constant timestep 'interval' and
probability of collision-coalescence determined by 'collcoalprob' */
{
  const double DELT(int2realtime(interval));

  const DoCoalescence coal{}; 

  const DoCollisions<Probability, DoCoalescence> colls(DELT,
                                                       collcoalprob,
                                                       coal);
  return ConstTstepMicrophysics(interval, colls);
}

#endif // COLLISIONCOALESCENCE_HPP