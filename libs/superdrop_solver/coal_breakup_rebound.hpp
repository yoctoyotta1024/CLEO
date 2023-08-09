// Author: Clara Bayley
// File: coal_breakup_rebound.hpp
/* Header file for class that enacts
collision events in which either coalescence, 
breakup, or rebound can occur.
CoalBreakupRebound struct satisfies
SDPairEnactX concept used in CollisionX struct */

#ifndef COAL_BREAKUP_REBOUND_HPP
#define COAL_BREAKUP_REBOUND_HPP

#include <functional>
#include <concepts>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./collisionxkernels.hpp"
#include "./collisionx.hpp"
#include "./coalescence.hpp"
#include "./breakup.hpp"

namespace dlc = dimless_constants;

class CoalBreakupRebound
/* class is method for coalescence / breakup between
two superdroplets. (Can be used in collisionsx struct
to enact collision-coalescence or collision-breakup
events in the superdroplet model) */
{
private:
  Coalescence coal;
  Breakup breakup;

  bool do_coalescence()
  /* returns true if coalescence shoulf occur based on the kinetic
  arguments in section 2.2 of Szakáll and Urbich 2018 (neglecting
  grazing angle considerations). */
  {

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
  CoalBreakupRebound(const double infrags) 
  : coal(Coalescence{}), breakup(Breakup(infrags)) {}

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
      if do_coalescence // (weber < x)
      {
        coal.coalesce_superdroplet_pair(drop1, drop2, gamma);
      }
      else // do breakup ie. weber > x
      {
        breakup.breakup_superdroplet_pair(drop1, drop2, gamma); 
      }
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