// Author: Clara Bayley
// File: breakup.hpp
/* Header file for class that enacts
collision-breakup events in
superdroplet model. Breakup struct
satisfies SDPairEnactX concept used in
CollisionX struct. Probability calculations
are contained in structures that satisfy the
requirements of the SDPairProbability
concept also used by CollisionX struct */

#ifndef BREAKUP_HPP
#define BREAKUP_HPP

#include "./superdrop.hpp"

class Breakup
/* class is method for enacting collisional-breakup given
two superdroplets. (Can be used in collisionsx
struct to enact collision-breakup events in SDM) */
{
private:
  void superdroplet_pair_breakup(Superdrop &drop1, Superdrop &drop2,
                                 const unsigned long long gamma) const
  /* enact collisional-breakup of droplets by changing multiplicity,
  radius and solute mass of each superdroplet in a pair. Method created
  by Author (no citation yet available) */
  {
    if (drop1.eps - gamma * drop2.eps == 0)
    {
      twin_superdroplet_breakup(drop1, drop2, gamma);
    }

    else if (drop1.eps - gamma * drop2.eps > 0)
    {
      different_superdroplet_breakup(drop1, drop2, gamma);
    }

    else
    {
      std::string errormsg = "something undefined occured "
                             "during colllision-breakup" +
                             std::to_string(drop1.eps) + " < " +
                             std::to_string(gamma * (drop2.eps));
      throw std::invalid_argument(errormsg);
    }
  }

  unsigned int breakup_gamma(const unsigned long long eps1,
                             const unsigned long long eps2,
                             const double prob,
                             const double phi) const
  /* calculates value of gamma factor in Monte Carlo
  collision-breakup, adapted from gamma for collision-
  coalescence in Shima et al. 2009. Here is is assumed
  maximally 1 breakup event can occur (gamma = 0 or 1)
  irrespective of if scaled probability, prob, is > 1 */
  {
    if (phi < (prob - floor(prob)))
    {
      return 1;
    }
    else // if phi >= (prob - floor(prob))
    {
      return 0;
    }
  }

public:
  void operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const
  /* this operator is used as an "adaptor" for using Breakup
  as a function in CollisionsX that satistfies the SDPairEnactX
  concept */
  {
    /* 1. calculate gamma factor for collision-breakup  */
    const unsigned int gamma = breakup_gamma(drop1.eps,
                                             drop2.eps,
                                             prob, phi);

    /* 2. enact collision-breakup on pair
    of superdroplets if gamma is not zero */
    if (gamma != 0)
    {
      superdroplet_pair_breakup(drop1, drop2, gamma);
    }
  }
};

template <SDPairProbability CollisionXProbability>
SdmProcess auto CollisionBreakupProcess(const int interval,
                                        const std::function<double(int)> int2time,
                                        const CollisionXProbability p)
{
  const double realtstep = int2time(interval);

  CollisionX<CollisionXProbability, Breakup>
      bus(realtstep, p, Breakup{});

  return ConstTstepProcess{interval, bus};
}

#endif // BREAKUP_HPP