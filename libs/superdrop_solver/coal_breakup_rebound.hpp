// Author: Clara Bayley
// File: coal_breakup_rebound.hpp
/* Header file for class that enacts
collision events in which either coalescence,
breakup, or rebound can occur.
CoalBreakupRebound struct satisfies
SDPairEnactX concept used in CollisionX struct */

#ifndef COAL_BREAKUP_REBOUND_HPP
#define COAL_BREAKUP_REBOUND_HPP

#include <algorithm>
#include <functional>
#include <concepts>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./collisionxkernels.hpp"
#include "./collisionx.hpp"
#include "./coalescence.hpp"
#include "./breakup.hpp"

namespace dlc = dimless_constants;

template <VelocityFormula TerminalVelocity>
class CoalBreakupRebound
/* class is method for coalescence / breakup between
two superdroplets. (Can be used in collisionsx struct
to enact collision-coalescence or collision-breakup
events in the superdroplet model) */
{
private:
  Coalescence coal;
  Breakup breakup;
  CollisionKinetics<TerminalVelocity> ck;

  unsigned int which_collisiontype(Superdrop &drop1,
                                   Superdrop &drop2) const
  /* based on the kinetic arguments in section 2.2 of
  Szakáll and Urbich 2018 (neglecting grazing angle considerations),
  returns 0 if rebound should occur, 1 if coalescence
  and 2 if breakup  */
  {
    const double cke(ck.collision_kinetic_energy(drop1, drop2));

    auto compare = [](const Superdrop &dropA, const Superdrop &dropB)
    {
      return dropA.radius < dropB.radius; // returns true if epsA < epsB
    };
    const auto smalldrop = std::min(drop1, drop2, compare); // drop with smaller drop.radius

    if (cke < ck.surfenergy(smalldrop))
    {
      return 0; // rebound
    }
    else if (cke < ck.coal_surfenergy(drop1, drop2)) // ie. Weber number < 1
    {
      return 1; // coalescence
    }
    else
    {
      return 2; // breakup
    }
  }

  unsigned long long collision_gamma(const unsigned long long eps1,
                                     const unsigned long long eps2,
                                     const double prob,
                                     const double phi) const
  /* calculates value of gamma factor in Monte Carlo
  collision as in Shima et al. 2009 */
  {
    return coal.coalescence_gamma(eps1, eps2, prob, phi);
  }

public:
  CoalBreakupRebound(TerminalVelocity tv, const double infrags)
      : coal(Coalescence{}), breakup(Breakup(infrags)), ck(tv) {}

  void operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const
  /* this operator is used as an "adaptor" for using
  CoalBreakupRebound as a function in CollisionsX
  that satistfies the SDPairEnactX concept */
  {
    /* 1. calculate gamma factor for collision-coalescence  */
    const unsigned long long gamma = collision_gamma(drop1.eps,
                                                     drop2.eps,
                                                     prob, phi);

    /* 2. enact collision-coalescence on pair
    of superdroplets if gamma is not zero */
    if (gamma != 0)
    {
      const unsigned int colltype(which_collisiontype(drop1, drop2));

      if (colltype == 1) // (weber < x)
      {
        coal.coalesce_superdroplet_pair(drop1, drop2, gamma);
      }

      else if (colltype == 2) // do breakup ie. weber > x
      {
        breakup.breakup_superdroplet_pair(drop1, drop2, gamma);
      }
      // else rebound, no change in superdroplets
    }
  }
};

template <SDPairProbability CollisionXProbability>
SdmProcess auto
CollisionCoalescenceProcess(const int interval,
                            const std::function<double(int)> int2time,
                            const CollisionXProbability p,
                            const double nfrags)
{
  const double realtstep = int2time(interval);

  CollisionX<CollisionXProbability, CoalBreakupRebound>
      coal(realtstep, p, CoalBreakupRebound(nfrags));

  return ConstTstepProcess{interval, coal};
};

#endif // COAL_BREAKUP_REBOUND_HPP