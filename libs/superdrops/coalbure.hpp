/*
 * ----- CLEO -----
 * File: coalbure.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 29th December 2023
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

struct CoalBuReFlag
{
  KOKKOS_FUNCTION
  unsigned int operator()(Superdrop &drop1,
                          Superdrop &drop2) const;
  /* function returns flag indicating rebound or
  coalescence or breakup. If flag = 1 -> coalescence.
  If flag = 2 -> breakup. Otherwise -> rebound.
  Flag decided based on the kinetic arguments in
  section 2.2 of Szakáll and Urbich 2018
  (neglecting grazing angle considerations) */
};

template <NFragments NFrags>
struct DoCoalBuRe
/* ie. DoCoalescenceBreakupRebound */
{
private:
  DoCoalescence coal;
  DoBreakup<NFrags> bu;
  CoalBuReFlag coalbure_flag;

  KOKKOS_FUNCTION
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

  KOKKOS_FUNCTION
  bool coalesce_breakup_or_rebound(const unsigned long long gamma,
                                   Superdrop &drop1,
                                   Superdrop &drop2) const;
  /*  function enacts rebound or coalescence or breakup
  depending on value of flag. If flag = 1 -> coalescence.
  If flag = 2 -> breakup. Otherwise -> rebound. */

public:
  DoCoalBuRe(const NFrags nfrags)
      : bu(nfrags), coalbure_flag(CoalBuReFlag{}) {}

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
  const auto DELT = double{int2realtime(interval)};

  const DoCoalBuRe coalbure(nfrags);
  const DoCollisions<Probability, DoCoalBuRe<NFrags>> colls(DELT,
                                                            collprob,
                                                            coalbure);

  return ConstTstepMicrophysics(interval, colls);
}

template <NFragments NFrags>
KOKKOS_FUNCTION bool
DoCoalBuRe<NFrags>::operator()(Superdrop &drop1, Superdrop &drop2,
                               const double prob, const double phi) const
/* this operator is used as an "adaptor" for
using DoCoalBuRe for collision - coalescence,
breakup or rebound as a function in DoCollisions
that satistfies the PairEnactX concept */
{
  /* 1. calculate gamma factor for collision  */
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();
  const auto gamma = collision_gamma(xi1, xi2, prob, phi);

  /* 2. enact collision between pair
  of superdroplets if gamma is not zero */
  if (gamma != 0)
  {
    return coalesce_breakup_or_rebound(gamma, drop1, drop2);
  }

  return 0;
}

template <NFragments NFrags>
KOKKOS_FUNCTION bool
DoCoalBuRe<NFrags>::
    coalesce_breakup_or_rebound(const unsigned long long gamma,
                                Superdrop &drop1,
                                Superdrop &drop2) const
/*  function enacts rebound or coalescence or breakup
depending on value of flag. If flag = 1 -> coalescence.
If flag = 2 -> breakup. Otherwise -> rebound. */
{
  const auto flag = coalbure_flag(drop1, drop2);

  bool is_null(0);
  switch (flag)
  {
  case 1: // coalescence
    is_null = coal.coalesce_superdroplet_pair(gamma, drop1, drop2);
    break;
  case 2: // breakup
    bu.breakup_superdroplet_pair(drop1, drop2);
    break;
  }

  return is_null;
}

KOKKOS_FUNCTION unsigned int
CoalBuReFlag::operator()(Superdrop &drop1,
                         Superdrop &drop2) const
/*  function returns flag indicating rebound or
coalescence or breakup. If flag = 1 -> coalescence.
If flag = 2 -> breakup. Otherwise -> rebound.
Flag decided based on the kinetic arguments in
section 2.2 of Szakáll and Urbich 2018
(neglecting grazing angle considerations) */
{
  const auto r1 = drop1.get_radius();
  const auto r2 = drop2.get_radius();
  const auto terminalv = SimmelTerminalVelocity{};

  const auto cke = collision_kinetic_energy(r1, r2,
                                            terminalv(drop1),
                                            terminalv(drop2));

  if (cke < surfenergy(Kokkos::fmin(r1, r2))) // cke < surface energy of small drop
  {
    return 0; // rebound
  }
  else if (cke < coal_surfenergy(r1, r2)) // weber number < 1 : coalescence
  {
    return 1; // coalescence
  }
  else // Weber > 1 : breakup
  {
    return 2; // breakup
  }
}

#endif // COALBURE_HPP