/*
 * ----- CLEO -----
 * File: breakup.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 23rd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality to enact collision- 
 * coalescence, breakup or rebound events
 * in SDM analagous to to Shima et al. 2009.
 * CoalBuRe struct satisfies PairEnactX
 * concept used in Collisions struct
 */

#ifndef COALBURE_HPP
#define COALBURE_HPP

#include <functional>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "./breakup.hpp"
#include "./coalescence.hpp"
#include "./collisions.hpp"
#include "./collisionkinetics.hpp"
#include "./microphysicalprocess.hpp"
#include "./superdrop.hpp"

template <NFragments NFrags>
struct DoCoalBuRe 
/* ie. DoCoalescenceBreakupRebound */ 
{
private:
  DoCoalescence coal;
  DoBreakup<NFrags> breakup;

public:
  DoCoalBuRe(const NFrags nfrags) : breakup(nfrags) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const;
  /* this operator is used as an "adaptor" for
  using DoCoalBuRe for collision - coalescence,
  breakup or rebound as a function in DoCollisions
  that satistfies the PairEnactX concept */
};

template <PairProbability Probability, NFragments NFrags>
inline MicrophysicalProcess auto
CoalBuRe(const unsigned int interval,
             const std::function<double(unsigned int)> int2realtime,
             const Probability collprob,
             const NFrags nfrags)
/* constructs Microphysical Process for collision-
coalscence, breakup or rebound of superdroplets with
a constant timestep 'interval' and probability
of collision determined by 'collprob' */
{
  const double DELT(int2realtime(interval));

  const DoCoalBuRe coalbure(nfrags);
  const DoCollisions<Probability, DoBreakup<NFrags>> colls(DELT,
                                                           collprob,
                                                           coalbure);

  return ConstTstepMicrophysics(interval, colls);
}

template <NFragments NFrags>
KOKKOS_INLINE_FUNCTION bool
DoCoalBuRe<NFrags>::operator()(Superdrop &drop1, Superdrop &drop2,
                               const double prob, const double phi) const
/* this operator is used as an "adaptor" for
using DoCoalBuRe for collision - coalescence,
breakup or rebound as a function in DoCollisions
that satistfies the PairEnactX concept */
{
  /* 1. calculate gamma factor for collision  */
  const unsigned long long xi1(drop1.get_xi());
  const unsigned long long xi2(drop2.get_xi());
  const unsigned long long gamma(collision_gamma(xi1, xi2, prob, phi));

  /* 2. enact collision between pair
  of superdroplets if gamma is not zero */
  if (gamma != 0)
  {
    coalesce_breakup_or_rebound(drop1, drop2, gamma);
  }
}

#endif // COALBURE_HPP