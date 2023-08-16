// Author: Clara Bayley
// File: coalescence.hpp
/* Header file for class that enacts
collision-coalescence events in
superdroplet model. Coalescence struct
satisfies SDinGBxPairEnactX concept used in
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
  void remove_empty_superdrop(SuperdropWithGbxindex &SDinGBx) const
  /* if multiplicity of drop = 0, ie. superdrop is empty,
  change it's SDinGBx to be value that indicates 
  superdrop is out of domain (ie. no longer exists) */
  {
    if (SDinGBx.superdrop.eps < 1) // ie. eps == 0
    {
      SDinGBx.sd_gbxindex = dlc::OUTOFDOMAIN;
    }
  }

  void twin_superdroplet_coalescence(SuperdropWithGbxindex &SDinGBx1,
                                     SuperdropWithGbxindex &SDinGBx2,
                                     const unsigned long long gamma) const
  /* if eps1 = gamma*eps2 coalescence makes twin SDs
  with same eps, r and solute mass. According to Shima et al. 2009
  Section 5.1.3. part (5) option (b)  */
  {
    Superdrop &sd1(SDinGBx1.superdrop);
    Superdrop &sd2(SDinGBx2.superdrop);

    const unsigned long long old_eps = sd2.eps; // = sd1.eps
    const unsigned long long new_eps = old_eps / 2;

    const double r1cubed(sd1.radius * sd1.radius * sd1.radius);
    const double r2cubed(sd2.radius * sd2.radius * sd2.radius);
    const double new_rcubed = r2cubed + gamma * r1cubed;
    const double new_r = std::pow(new_rcubed, (1.0 / 3.0));

    const double new_m_sol = sd2.m_sol + gamma * sd1.m_sol;

    sd1.eps = new_eps;
    sd2.eps = old_eps - new_eps;

    sd1.radius = new_r;
    sd2.radius = new_r;

    sd1.m_sol = new_m_sol;
    sd2.m_sol = new_m_sol;
  
    remove_empty_superdrop(SDinGBx1); // because if eps1 = eps2 = 1 before coalesence, then eps1=0 now
  }

  void different_superdroplet_coalescence(Superdrop &sd1,
                                          Superdrop &sd2,
                                          const unsigned long long gamma) const
  /* if eps1 > gamma*eps2 coalescence grows sd2 radius and mass
  via decreasing multiplicity of sd1. According to
  Shima et al. 2009 Section 5.1.3. part (5) option (a)  */
  {
    sd1.eps = sd1.eps - gamma * sd2.eps;

    const double r1cubed(sd1.radius * sd1.radius * sd1.radius);
    const double r2cubed(sd2.radius * sd2.radius * sd2.radius);
    const double new_rcubed = r2cubed + gamma * r1cubed;

    sd2.radius = std::pow(new_rcubed, (1.0 / 3.0));
    sd2.m_sol = sd2.m_sol + gamma * sd1.m_sol;
  }

public:
  void operator()(SuperdropWithGbxindex &SDinGBx1,
                  SuperdropWithGbxindex &SDinGBx2,
                  const double prob, const double phi) const
  /* this operator is used as an "adaptor" for using Coalescence
  as a function in CollisionsX that satistfies the SDinGBxPairEnactX
  concept */
  {
    const unsigned long long eps1(SDinGBx1.superdrop.eps);
    const unsigned long long eps2(SDinGBx2.superdrop.eps);

    /* 1. calculate gamma factor for collision-coalescence  */
    const unsigned long long gamma(coalescence_gamma(eps1, eps2,
                                                     prob, phi));

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
    const unsigned long long eps1(SDinGBx1.superdrop.eps);
    const unsigned long long eps2(SDinGBx2.superdrop.eps);

    if (eps1 - gamma * eps2 > 0)
    {
      different_superdroplet_coalescence(SDinGBx1.superdrop,
                                         SDinGBx2.superdrop,
                                         gamma);
    }

    else if (eps1 - gamma * eps2 == 0)
    {
      twin_superdroplet_coalescence(SDinGBx1, SDinGBx2, gamma);
    }

    else
    {
      std::string errormsg("something undefined occured "
                           "during colllision-coalescence");
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
