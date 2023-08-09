// Author: Clara Bayley
// File: coal_breakup.hpp
/* Header file for class that enacts
collision events in which either coalescence 
or breakup can occur. CoalBreakup struct
satisfies SDPairEnactX concept used in
CollisionX struct */

#ifndef COAL_BREAKUP_HPP
#define COAL_BREAKUP_HPP

#include <functional>
#include <concepts>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./collisionxkernels.hpp"
#include "./collisionx.hpp"

namespace dlc = dimless_constants;

class CoalBreakup
/* class is method for coalescence / breakup between
two superdroplets. (Can be used in collisionsx struct
to enact collision-coalescence or collision-breakup
events in the superdroplet model) */
{
private:




  unsigned long long collision_gamma(const unsigned long long eps1,
                                       const unsigned long long eps2,
                                       const double prob,
                                       const double phi) const
  /* calculates value of gamma factor in Monte Carlo
  collision-coalescence as in Shima et al. 2009 */
  {
    unsigned long long gamma = floor(prob); // if phi >= (prob - floor(prob))
    if (phi < (prob - gamma))
    {
      ++gamma;
    }

    const unsigned long long maxgamma(eps1 / eps2); // same as floor() for positive ints

    return std::min(gamma, maxgamma);
  }

public:
  void operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const
  /* this operator is used as an "adaptor" for using CoalBreakup
  as a function in CollisionsX that satistfies the SDPairEnactX
  concept */
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
        coalesce_superdroplet_pair(drop1, drop2, gamma);
      }
      else // do breakup ie. weber > x
      {
        breakup_superdroplet_pair(drop1, drop2, gamma); 
      }
    } 
  }
};

template <SDPairProbability CollisionXProbability>
SdmProcess auto
CollisionCoalescenceProcess(const int interval,
                            const std::function<double(int)> int2time,
                            const CollisionXProbability p)
{
  const double realtstep = int2time(interval);

  CollisionX<CollisionXProbability, CoalBreakup>
      coal(realtstep, p, CoalBreakup{});

  return ConstTstepProcess{interval, coal};
};

#endif // COAL_BREAKUP_HPP