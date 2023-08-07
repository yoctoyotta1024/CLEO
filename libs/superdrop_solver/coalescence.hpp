// Author: Clara Bayley
// File: coalescence.hpp
/* Header file for class that enacts
collision-coalescence events in
superdroplet model. Coalescence struct
satisfies SDPairEnactX concept used in
CollisionX struct. Probability calculations
are contained in structures that satisfy the
requirements of the SDPairProbability
concept also used by CollisionX struct */

#ifndef COALESCENCE_HPP
#define COALESCENCE_HPP

#include <string>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <concepts>
#include <cmath>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./collisionxkernels.hpp"
#include "./collisionx.hpp"

namespace dlc = dimless_constants;

class Coalescence
/* class is method for coalescence between
two superdroplets. (Can be used in collisionsx struct
to enact collision-coalescence events in SDM) */
{
private:
  void coalesce_superdroplet_pair(Superdrop &drop1, Superdrop &drop2,
                                  const unsigned long long gamma) const
  /* coalesce pair of superdroplets by changing multiplicity,
  radius and solute mass of each superdroplet in pair
  according to Shima et al. 2009 Section 5.1.3. part (5) */
  {
    if (drop1.eps - gamma * drop2.eps == 0)
    {
      twin_superdroplet_coalescence(drop1, drop2, gamma);
    }

    else if (drop1.eps - gamma * drop2.eps > 0)
    {
      different_superdroplet_coalescence(drop1, drop2, gamma);
    }

    else
    {
      std::string errormsg = "something undefined occured "
                             "during colllision-coalescence" +
                             std::to_string(drop1.eps) + " < " +
                             std::to_string(gamma * (drop2.eps));
      throw std::invalid_argument(errormsg);
    }
  }

  void twin_superdroplet_coalescence(Superdrop &drop1,
                                     Superdrop &drop2,
                                     const unsigned long long gamma) const
  /* if eps1 = gamma*eps2 coalescence makes twin SDs
  with same eps, r and solute mass. According to Shima et al. 2009
  Section 5.1.3. part (5) option (b)  */
  {
    const unsigned long long new_eps = drop2.eps / 2;
    const double new_m_sol = drop2.m_sol + gamma * drop1.m_sol;
    const double new_rcubed = pow(drop2.radius, 3.0) + gamma * pow(drop1.radius, 3.0);
    const double new_r = pow(new_rcubed, (1.0 / 3.0));

    drop1.eps = new_eps;
    drop2.eps = drop2.eps - new_eps;

    drop1.radius = new_r;
    drop2.radius = new_r;

    drop1.m_sol = new_m_sol;
    drop2.m_sol = new_m_sol;
  }

  void different_superdroplet_coalescence(Superdrop &drop1,
                                          Superdrop &drop2,
                                          const unsigned long long gamma) const
  /* if eps1 > gamma*eps2 coalescence grows drop2 radius and mass
  via decreasing multiplicity of drop1. According to
  Shima et al. 2009 Section 5.1.3. part (5) option (a)  */
  {
    drop1.eps = drop1.eps - gamma * drop2.eps;

    const double new_rcubed = pow(drop2.radius, 3.0) + gamma * pow(drop1.radius, 3.0);
    drop2.radius = pow(new_rcubed, (1.0 / 3.0));
    drop2.m_sol = drop2.m_sol + gamma * drop1.m_sol;
  }

  unsigned long long coalescence_gamma(const unsigned long long eps1,
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
  /* this operator is used as an "adaptor" for using Coalescence
  as a function in CollisionsX that satistfies the SDPairEnactX
  concept */
  {
    /* 1. calculate gamma factor for collision-coalescence  */
    const unsigned long long gamma = coalescence_gamma(drop1.eps,
                                                       drop2.eps,
                                                       prob, phi);

    /* 2. enact collision-coalescence on pair
    of superdroplets if gamma is not zero */
    if (gamma != 0)
    {
      coalesce_superdroplet_pair(drop1, drop2, gamma);
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

  CollisionX<CollisionXProbability, Coalescence>
      coal(realtstep, p, Coalescence{});

  return ConstTstepProcess{interval, coal};
};

#endif // COALESCENCE_HPP
