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
#include "./terminalvelocity.hpp"
#include "./collisionxkernels.hpp"
#include "./collisionx.hpp"
#include "./coalescence.hpp"
#include "./breakup.hpp"

namespace dlc = dimless_constants;

template <VelocityFormula TerminalVelocity>
class CoalBreakupRebound
/* class is method for coalescence / breakup / rebound
between two superdroplets. (Can be used in collisionsx
struct to enact collision-coalescence or collision-breakup
or collision-rebound events in the superdroplet model) */
{
private:
  Coalescence coal;
  Breakup breakup;
  CollisionKinetics<TerminalVelocity> ck;

  void coalesce_or_breakup(Superdrop &sd1,
                           Superdrop &sd2,
                           const unsigned long long gamma) const
  /* based on the kinetic arguments in section 2.2 of
  Szakáll and Urbich 2018 (neglecting grazing angle considerations),
  function enacts coalescence or breakup */
  {
    const double cke(ck.collision_kinetic_energy(sd1, sd2));

    if (cke < ck.coal_surfenergy(sd1, sd2)) // Weber number < 1 : coalescence
    {
      coal.coalesce_superdroplet_pair(sd1, sd2, gamma);
    }
    else // Weber > 1 : breakup
    {
      breakup.breakup_superdroplet_pair(sd1, sd2);
    }
  }

  void coalesce_breakup_or_rebound(Superdrop &sd1,
                                   Superdrop &sd2,
                                   const unsigned long long gamma) const
  /* based on the kinetic arguments in section 2.2 of
  Szakáll and Urbich 2018 (neglecting grazing angle considerations),
  function enacts rebound or coalescence/breakup */
  {
    auto compare = [](const Superdrop &dropA, const Superdrop &dropB)
    {
      return dropA.radius < dropB.radius; // returns true if epsA < epsB
    };
    const auto smalldrop(std::min(sd1, sd2, compare)); // drop with smaller drop.radius

    const double cke(ck.collision_kinetic_energy(sd1, sd2));

    if (cke >= ck.surfenergy(smalldrop)) // ie. if not rebound
    {
      coalesce_or_breakup(sd1, sd2, gamma);
    } 
  }

  unsigned long long collision_gamma(const unsigned long long eps1,
                                     const unsigned long long eps2,
                                     const double probcoll,
                                     const double phi) const
  /* calculates value of gamma factor in Monte Carlo
  collision as in Shima et al. 2009 given probability of
  collision.
  *note* argument is NOT probability of collision-coalescence! */
  {
    return coal.coalescence_gamma(eps1, eps2, probcoll, phi);
  }

public:
  CoalBreakupRebound(TerminalVelocity tv, const double infrags)
      : coal(Coalescence{}), breakup(Breakup(infrags)), ck(tv) {}

  void operator()(SuperdropWithGbxindex &SDinGBx1,
                  SuperdropWithGbxindex &SDinGBx2,
                  const double probcoll, const double phi) const
  /* this operator is used as an "adaptor" for using
  CoalBreakupRebound as a function in CollisionsX
  that satistfies the SDPairEnactX concept.
  *note* operator uses probcoll, probability of collision,
  NOT probability of collision-coalescence! */
  {
    Superdrop &sd1(SDinGBx1.superdrop);
    Superdrop &sd2(SDinGBx2.superdrop);

    /* 1. calculate gamma factor for collision  */
    const unsigned long long gamma(collision_gamma(sd1.eps,
                                                   sd2.eps,
                                                   probcoll, phi));

    /* 2. enact collision between pair
    of superdroplets if gamma is not zero */
    if (gamma != 0)
    {
      coalesce_breakup_or_rebound(sd1, sd2, gamma);    // include rebound
      // coalesce_or_breakup(sd1, sd2, gamma);         // no rebound
    }
  }
};

template <SDPairProbability CollisionXProbability,
          VelocityFormula TerminalVelocity>
SdmProcess auto
CollisionAllProcess(const int interval,
                    const std::function<double(int)> int2time,
                    const CollisionXProbability p,
                    const TerminalVelocity tv,
                    const double nfrags)
/* SDM process for collisions of superdroplets 
followed by coalescence, breakup or rebound */
{
  const double realtstep = int2time(interval);

  CollisionX<CollisionXProbability,
             CoalBreakupRebound<TerminalVelocity>>
      collall(realtstep, p, CoalBreakupRebound(tv, nfrags));

  return ConstTstepProcess{interval, collall};
};

#endif // COAL_BREAKUP_REBOUND_HPP