/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: coalbure.hpp
 * Project: collisionprobs
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 11th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality to enact collision-
 * coalescence, breakup or rebound events
 * in SDM analagous to to Shima et al. 2009.
 * CoalBuRe struct satisfies PairEnactX
 * concept used in Collisions struct
 */


#ifndef LIBS_SUPERDROPS_COLLISIONPROBS_COALBURE_HPP_
#define LIBS_SUPERDROPS_COLLISIONPROBS_COALBURE_HPP_

#include <concepts>
#include <functional>

#include <Kokkos_Core.hpp>

#include "./breakup.hpp"
#include "./breakup_nfrags.hpp"
#include "./coalbure_flag.hpp"
#include "./coalescence.hpp"
#include "./collisionkinetics.hpp"
#include "./collisions.hpp"
#include "../microphysicalprocess.hpp"
#include "../superdrop.hpp"

/* ie. DoCoalescenceBreakupRebound */
template <NFragments NFrags, CoalBuReFlag Flag>
struct DoCoalBuRe {
 private:
  DoCoalescence coal;
  DoBreakup<NFrags> bu;
  Flag coalbure_flag;

  /* calculates value of gamma factor in Monte Carlo
  collision as in Shima et al. 2009 given probability of
  collision. Note: probability is probability of
  collision *NOT* collision-coalescence! */
  KOKKOS_FUNCTION
  uint64_t collision_gamma(const uint64_t xi1, const uint64_t xi2, const double prob,
                           const double phi) const {
    return coal.coalescence_gamma(xi1, xi2, prob, phi);
  }

  /*  function enacts rebound or coalescence or breakup
  depending on value of flag. If flag = 1 -> coalescence.
  If flag = 2 -> breakup. Otherwise -> rebound. */
  KOKKOS_FUNCTION
  bool coalesce_breakup_or_rebound(const uint64_t gamma, const double phi, Superdrop &drop1,
                                   Superdrop &drop2) const;

 public:
  DoCoalBuRe(const NFrags nfrags, const Flag flag) : bu(nfrags), coalbure_flag(flag) {}

  /* this operator is used as an "adaptor" for
  using DoCoalBuRe for collision - coalescence,
  breakup or rebound as a function in DoCollisions
  that satistfies the PairEnactX concept */
  KOKKOS_INLINE_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2, const double prob, const double phi) const;
};

/* constructs Microphysical Process for collision-
coalscence, breakup or rebound of superdroplets with
a constant timestep 'interval' and probability
of collision determined by 'collprob' */
template <PairProbability Probability, NFragments NFrags, CoalBuReFlag Flag>
inline MicrophysicalProcess auto CoalBuRe(const unsigned int interval,
                                          const std::function<double(unsigned int)> int2realtime,
                                          const Probability collprob, const NFrags nfrags,
                                          const Flag coalbure_flag) {
  const auto DELT = double{int2realtime(interval)};

  const DoCoalBuRe<NFrags, Flag> coalbure(nfrags, coalbure_flag);
  const DoCollisions<Probability, DoCoalBuRe<NFrags, Flag>> colls(DELT, collprob, coalbure);

  return ConstTstepMicrophysics(interval, colls);
}

/* this operator is used as an "adaptor" for
using DoCoalBuRe for collision - coalescence,
breakup or rebound as a function in DoCollisions
that satistfies the PairEnactX concept */
template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool DoCoalBuRe<NFrags, Flag>::operator()(Superdrop &drop1, Superdrop &drop2,
                                                          const double prob,
                                                          const double phi) const {
  /* 1. calculate gamma factor for collision  */
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();
  const auto gamma = collision_gamma(xi1, xi2, prob, phi);

  /* 2. enact collision between pair
  of superdroplets if gamma is not zero */
  if (gamma != 0) {
    return coalesce_breakup_or_rebound(gamma, phi, drop1, drop2);
  }

  return 0;
}

/*  function enacts rebound or coalescence or breakup
depending on value of flag. If flag = 1 -> coalescence.
If flag = 2 -> breakup. Otherwise -> rebound. */
template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool DoCoalBuRe<NFrags, Flag>::coalesce_breakup_or_rebound(const uint64_t gamma,
                                                                           const double phi,
                                                                           Superdrop &drop1,
                                                                           Superdrop &drop2) const {
  const auto flag = coalbure_flag(phi, drop1, drop2);

  bool is_null(0);
  switch (flag) {
    case 1:  // coalescence
      is_null = coal.coalesce_superdroplet_pair(gamma, drop1, drop2);
      break;
    case 2:  // breakup
      bu.breakup_superdroplet_pair(drop1, drop2);
      break;
  }

  return is_null;
}

#endif  // LIBS_SUPERDROPS_COLLISIONPROBS_COALBURE_HPP_
