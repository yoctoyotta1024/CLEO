// Author: Clara Bayley
// File: coalescence.hpp
/* Header file for class that enacts
collision-coalescence events in
superdroplet model. Coalescence struct
satisfies SDPairEnactX concept used in
CollisionX struct */

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
  void remove_null_superdrops(Superdrop &drop2)
  {
    sd.gbxindex = limits; //TODO
  }

  void twin_superdroplet_coalescence(SuperdropWithGbxindex &SDinGBx1,
                                     SuperdropWithGbxindex &SDinGBx2,
                                     const unsigned long long gamma) const
  /* if eps1 = gamma*eps2 coalescence makes twin SDs
  with same eps, r and solute mass. According to Shima et al. 2009
  Section 5.1.3. part (5) option (b)  */
  {
    Superdrop &d1(SDinGBx1.superdrop);
    Superdrop &d2(SDinGBx2.superdrop);

    const unsigned long long new_eps(d2.eps / 2); // = d1.eps /2 
    d1.eps = new_eps
    d2.eps = d2.eps - new_eps;

    const double r2cubed(d2.radius * d2.radius * d2.radius);
    const double r1cubed(d1.radius * d1.radius * d1.radius);
    const double new_rcubed(r2cubed + gamma * r1cubed);
    const double new_r(std::pow(new_rcubed, (1.0 / 3.0)));
    d1.radius = new_r;
    d2.radius = new_r;
    
    const double new_m_sol = d2.m_sol + gamma * d1.m_sol;
    d1.m_sol = new_m_sol;
    d2.m_sol = new_m_sol;

    remove_null_superdrops(drop2);
  }

  void different_superdroplet_coalescence(Superdrop &d1,
                                          Superdrop &d2,
                                          const unsigned long long gamma) const
  /* if eps1 > gamma*eps2 coalescence grows d2 radius and mass
  via decreasing multiplicity of d1. According to
  Shima et al. 2009 Section 5.1.3. part (5) option (a)  */
  {
    d1.eps = d1.eps - gamma * d2.eps;

    const double r2cubed(d2.radius * d2.radius * d2.radius);
    const double r1cubed(d1.radius * d1.radius * d1.radius);
    const double new_rcubed(r2cubed + gamma * r1cubed);
    d2.radius = std::pow(new_rcubed, (1.0 / 3.0));
 
    d2.m_sol = d2.m_sol + gamma * d1.m_sol;
  }

public:
  void operator()(SuperdropWithGbxindex &SDinGBx1,
                  SuperdropWithGbxindex &SDinGBx2,
                  const double prob, const double phi) const
  /* this operator is used as an "adaptor" for using Coalescence
  as a function in CollisionsX that satistfies the SDPairEnactX
  concept */
  {
    /* 1. calculate gamma factor for collision-coalescence  */
    const unsigned long long eps1(SDinGBx1.superdrop.eps);
    const unsigned long long eps2(SDinGBx2.superdrop.eps);
    const unsigned long long gamma = coalescence_gamma(eps1, eps2,
                                                       prob, phi);

    /* 2. enact collision-coalescence on pair
    of superdroplets if gamma is not zero */
    if (gamma != 0)
    {
      coalesce_superdroplet_pair(SDinGBx1, SDinGBx2, gamma);
    }
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

  void coalesce_superdroplet_pair(SuperdropWithGbxindex &SDinGBx1,
                                  SuperdropWithGbxindex &SDinGBx2,
                                  const unsigned long long gamma) const
  /* coalesce pair of superdroplets by changing multiplicity,
  radius and solute mass of each superdroplet in pair
  according to Shima et al. 2009 Section 5.1.3. part (5) */
  {
    if (drop1.eps - gamma * drop2.eps > 0)
    {
      different_superdroplet_coalescence(SDinGBx1.superdrop,
                                         SDinGBx2.superdrop, gamma);
    }

    else if (drop1.eps - gamma * drop2.eps == 0)
    {
      twin_superdroplet_coalescence(SDinGBx1, SDinGBx2, gamma);
    }

    else
    {
      std::string errormsg = "something undefined occured "
                             "during colllision-coalescence";
      throw std::invalid_argument(errormsg);
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
