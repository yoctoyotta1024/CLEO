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

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./hydrodynamicprob.hpp"
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

struct GolovinCollCoalProb
/* Probability of collision-coalescence of
a pair of droplets according to Golovin 1963
(see e.g. Shima et al. 2009) */
{
  const double prob_jk_const;

  GolovinCollCoalProb()
      : prob_jk_const(1.5e3 * (pow(dlc::R0, 3.0))) {}

  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const
  /* returns probability that a pair of droplets coalesces
  according to Golovin's (sum of volumes) coalescence kernel.
  Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol where
  K(drop1, drop2) := C(drop1, drop2) * |v1−v2|, (see Shima 2009 eqn 3),
  and K(drop1, drop2) is Golovin 1963 (coalescence) kernel */
  {
    const double DELT_DELVOL = DELT / VOLUME;                                   // time interval / volume for which collision probability is calculated [s/m^3]
    const double golovins_kernel = prob_jk_const * (drop1.vol() + drop2.vol()); // Golovin 1963 coalescence kernel

    const double prob_jk = golovins_kernel * DELT_DELVOL;

    return prob_jk;
  }
};

struct LongKernelEfficiency
/* Collision-Coalescence Efficiency factor in Long's
  Hydrodynamic kernel according to Simmel et al. 2002.
  eff = collision-coalescence efficiency E(R,r) where R>r.
  eff = colleff(R,r) * coaleff(R,r) (see eqn 12 Simmel et al. 2002).
  Here it's assumed that coaleff(R,r) = 1, which also means that
  for collisions where R > rlim, eff(R,r) = colleff(R,r) = 1. */
{
  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2) const
  {
    constexpr double coaleff = 1.0;
    constexpr double rlim = 5e-5 / dlc::R0;          // 50 micron limit to determine collision-coalescence efficiency (eff)
    constexpr double colleff_lim = 0.001;            // minimum efficiency if larger droplet's radius < rlim
    constexpr double A1 = 4.5e4 * dlc::R0 * dlc::R0; // constants in efficiency calc if larger droplet's radius < rlim
    constexpr double A2 = 3e-4 / dlc::R0;

    const double bigr = std::max(drop1.radius, drop2.radius);
    const double smallr = std::min(drop1.radius, drop2.radius);

    /* calculate collision-coalescence efficiency, eff = colleff * coaleff */
    double colleff(1.0);
    if (bigr < rlim)
    {
      const double colleff_calc = A1 * pow(bigr, 2.0) * (1 - A2 / smallr);
      colleff = std::max(colleff_calc, colleff_lim); // colleff >= colleff_lim
    }

    const double eff = colleff * coaleff;

    return eff;
  }
};

HydrodynamicProb<LongKernelEfficiency, SimmelTerminalVelocity>
LongCollCoalProb()
/* returns the probability of collision-coalescence
using Simmel et al. 2002's formulation of
Long's Hydrodynamic Kernel */
{
  return HydrodynamicProb(LongKernelEfficiency{},
                          SimmelTerminalVelocity{});
}

template <SDPairProbability CollisionXProbability>
SdmProcess auto CollisionCoalescenceProcess(const int interval,
                                            const std::function<double(int)> int2time,
                                            const CollisionXProbability p)
{
  const double realtstep = int2time(interval);

  CollisionX<CollisionXProbability, Coalescence>
      coals(realtstep, p, Coalescence{});

  return ConstTstepProcess{interval, coals};
}

#endif // COALESCENCE_HPP