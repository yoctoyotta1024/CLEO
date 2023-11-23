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

  KOKKOS_INLINE_FUNCTION
  unsigned long long collision_gamma(const unsigned long long xi1,
                                     const unsigned long long xi2,
                                     const double prob,
                                     const double phi) const
  /* calculates value of gamma factor in Monte Carlo
  collision as in Shima et al. 2009 given probability of
  collision. Note: probability is probability of
  collision *NOT* collision-coalescence! */
  {
    return coal.coalescence_gamma(xi1, xi2, prob, phi);
  }

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
    return coalesce_breakup_or_rebound(drop1, drop2, gamma);
  }

  return 0;
}

template <NFragments NFrags>
KOKKOS_INLINE_FUNCTION bool
DoCoalBuRe<NFrags>::coalesce_breakup_or_rebound(Superdrop &drop1,
                                                Superdrop &drop2,
                                                const unsigned long long gamma) const
/* based on the kinetic arguments in section 2.2 of
Szakáll and Urbich 2018 (neglecting grazing angle considerations),
function enacts rebound or coalescence/breakup */
{
  const Superdrop &sd1(SDinGBx1.superdrop);
  const Superdrop &sd2(SDinGBx2.superdrop);

  auto compare = [](const Superdrop &dropA, const Superdrop &dropB)
  {
    return dropA.radius < dropB.radius; // returns true if epsA < epsB
  };
  const auto smalldrop(std::min(sd1, sd2, compare)); // drop with smaller drop.radius

  const double cke(ck.collision_kinetic_energy(sd1, sd2));

  if (cke >= ck.surfenergy(smalldrop)) // ie. if not rebound
  {
    coalesce_or_breakup(SDinGBx1, SDinGBx2, gamma);
  } 
}

#endif // COALBURE_HPP