// Author: Clara Bayley
// File: breakup.hpp
/* Header file for class that enacts
collision-breakup events in
superdroplet model. Breakup struct
satisfies SDPairEnactX concept used in
CollisionX struct */

#ifndef BREAKUP_HPP
#define BREAKUP_HPP

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cmath>

#include "./superdrop.hpp"

class Breakup
/* class is method for enacting collisional-breakup given
two superdroplets. (Can be used in collisionsx
struct to enact collision-breakup events in SDM) */
{
private:
  double nfrags; // expected number of fragments from 1 real droplet collision (Omega in my notes)

  void superdroplet_pair_breakup(Superdrop &drop1, Superdrop &drop2) const
  /* enact collisional-breakup of droplets by changing multiplicity,
  radius and solute mass of each superdroplet in a pair. Method created
  by Author (no citation yet available). Note implicit assumption that
  gamma factor = 1. */
  {
    if (drop1.eps == drop2.eps)
    {
      twin_superdroplet_breakup(drop1, drop2);
    }

    else
    {
      different_superdroplet_breakup(drop1, drop2);
    }
  }

  void twin_superdroplet_breakup(Superdrop &drop1,
                                 Superdrop &drop2) const
  /* if eps1 = gamma*eps2 breakup of same multiplicity SDs
  produces (non-identical) twin SDs. Similar to
  Shima et al. 2009 Section 5.1.3. part (5) option (b).
  Note implicit assumption that gamma factor = 1. */
  {
    const unsigned long long old_eps = drop2.eps; // = drop1.eps
    const unsigned long long new_eps = std::round(nfrags * old_eps) / 2;
    const double sumr3 = std::pow(drop1.radius, 3.0) +
                             std::pow(drop2.radius, 3.0);
    const double new_r = std::pow(old_eps / new_eps * sumr3, (1.0/3.0));
    const double new_m_sol = old_eps * (drop1.m_sol + drop2.m_sol) / new_eps;
 
    drop1.eps = new_eps;
    drop2.eps = old_eps - new_eps;

    drop1.radius = new_r;
    drop2.radius = new_r;

    drop1.m_sol = new_m_sol;
    drop2.m_sol = new_m_sol;    
  }

  void different_superdroplet_breakup(Superdrop &drop1,
                                      Superdrop &drop2) const
  /* if eps1 > gamma*eps2 breakup alters drop2 radius and mass
  via decreasing multiplicity of drop1. Similar to
  Shima et al. 2009 Section 5.1.3. part (5) option (a).
  Note implicit assumption that gamma factor = 1. */
  {
    drop1.eps = drop1.eps - drop2.eps;

    const unsigned long long old_eps = drop2.eps;
    const unsigned long long new_eps = std::round(nfrags * old_eps);
    const double sumr3 = std::pow(drop1.radius, 3.0) +
                             std::pow(drop2.radius, 3.0);
    
    drop2.eps = new_eps;
    drop2.radius = std::pow(sumr3 * old_eps / new_eps, (1.0/3.0)); // NOTE: implicit casting of eps from unsigned long long to double here
    drop2.m_sol = old_eps * (drop1.m_sol + drop2.m_sol) / new_eps; 
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
  Breakup(const double infrags) : nfrags(std::max(infrags, 1.0))
  /* nfrags = expected number of fragments from one real
  droplet collision-breakup event (Omega in my notes). Conservative
  estimate is that nfrags >=1.0 to ensure largest possible 
  fragment is <= mass1 + mass2 (sum of original drop masses) */
  {
    if (infrags < 0.0)
    {
      const std::string err("attempted to initialise breakup"
                            "with invalid value for nfrags");
      throw std::invalid_argument(err);
    }
  }

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
      superdroplet_pair_breakup(drop1, drop2);
    }
  }
};

template <SDPairProbability CollisionXProbability>
SdmProcess auto
CollisionBreakupProcess(const int interval,
                        const std::function<double(int)> int2time,
                        const CollisionXProbability p,
                        const double nfrags)
{
  const double realtstep = int2time(interval);

  CollisionX<CollisionXProbability, Breakup>
      bu(realtstep, p, Breakup(nfrags));

  return ConstTstepProcess{interval, bu};
}

#endif // BREAKUP_HPP