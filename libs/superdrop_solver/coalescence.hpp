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
/* Collision-Coalescence Efficiency factor, eff, in Long's
  Hydrodynamic kernel according to Simmel et al. 2002.
  eff = collision-coalescence efficiency E(R,r) where R>r.
  eff = colleff(R,r) * coaleff(R,r) (see eqn 12 and 13
  Simmel et al. 2002). Here it's usually assumed that
  coaleff(R,r) = 1, ie. eff = colleff. (This would also mean
  that for collisions where R > rlim, eff(R,r) = colleff(R,r) = 1) */
{
  const double coaleff;

  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2) const
  {
    constexpr double rlim = 5e-5 / dlc::R0;          // 50 micron limit to determine collision-coalescence efficiency (eff)
    constexpr double colleff_lim = 0.001;            // minimum efficiency if larger droplet's radius < rlim
    constexpr double A1 = 4.5e4 * dlc::R0 * dlc::R0; // constants in efficiency calc if larger droplet's radius < rlim
    constexpr double A2 = 3e-4 / dlc::R0;

    const auto [smallr, bigr] = std::minmax(drop1.radius, drop2.radius);

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

template <VelocityFormula TerminalVelocity>
struct LowListKernelEfficiency
/* Collision-Coalescence Efficiency factor, eff, for the
  Hydrodynamic kernel. eff = colleff(R,r) * coaleff(R,r) where:
  - colleff is Long's collision efficiency as seen in equation 13
  of Simmel et al. 2002
  - coaleff is from equation (4.5) and (4.6) Low and List 1982(a) */
{
  TerminalVelocity terminalv;
  LongKernelEfficiency colleff{1.0};

  double kinetic_energy(const Superdrop &drop1,
                        const Superdrop &drop2)
  /* returns cke/ pi, where cke = collision kinetic energy
  as formulated in Low and List 1982(a) eqn 3.1 */
  {
    const double r1(drop1.radius);
    const double r2(drop2.radius);
    const double ratio = std::pow(r1, 3.0) / (1 + std::pow(r1/r2, 3.0));
    
    const double vdiff = (terminalv(drop1) - terminalv(drop2)) * dlc::W0; // [m/s]

    const double cke_pi = DC::RHO_L / 12.0 * ratio * vdiff * vdiff;

    return cke_pi;
  }

  double total_surfenergy(const Superdrop &drop1,
                          const Superdrop &drop2)
  /* returns total surface energy of drops divided by pi
  as in equation 4.2 of Low and List 1982 */
  {
    constexpr double sigma = 7.28e-2;  // [J/m^-2]
    constexpr double surfconst = 4.0 * sigma; // [J/m^-2]
    const double r1(drop1.radius) * dlc::R0; // [m]
    const double r2(drop2.radius) * dlc::R0; // [m] 

    const double r2sum = (r1 * r1 + r2 * r2);
    const double tot_surfe_pi =  surfconst * r2sum;
    
    return tot_surfe_pi // total surface energy / pi
  }

  double equivalent_surfenergy(const Superdrop &drop1,
                               const Superdrop &drop2)
  /* returns surface energy of single spherical equivalent of
  drops divided by pi as in equation 4.3 of Low and List 1982 */
  {
    constexpr double sigma = 7.28e-2;  // [J/m^-2]
    constexpr double surfconst = 4.0 * sigma; // [J/m^-2]
    const double r1(drop1.radius) * dlc::R0; // [m]
    const double r2(drop2.radius) * dlc::R0; // [m]  

    const double r3sum = std::pow(r1, 3.0) + std::pow(r2, 3.0);
    const double equiv_surfe_pi = surfconst * std::pow(r3sum, 2.0/3.0); 
    
    return equiv_surfe_pi // spherical equivalent surface energy / pi
  }

  double exponential(const double etot_pi, const double surf_c_pi)
  /* calc exponential in eqn 4.5 Low and List 1982(a) given
  total collision energy, etot/pi [J] and equivalent surface
  energy, surf_c / pi [J] */
  {
    constexpr double bconst = -2.62e6; // [J^-2]
    constexpr double sigma = 7.28e-2;  // [J/m^-2]

    const double exponent = bconst * sigma * std::numbers::pi *
                              etot_pi * etot_pi / surf_c_pi;

    return std::exp(exponent);
  }

  double sizeratio_factor(const double r1, const double r2)
  /* calc factor that takes into accoutn size ratio fo droplets in
  eqn 4.5 Low and List 1982(a). */
  {
    const double alpha(1 + std::min(r1, r2) / std::max(r1, r2)); // alpha = 1 + Ds/Dl 
    
    return 1.0 / (alpha * alpha); // alpha^(-2)
  }

  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2) const
  {
    constexpr double aconst = 0.778;
    constexpr double energylim = 5e-6 / std::numbers::pi ; // etot limit / pi [J]

    const double surf_t_pi = total_surfenergy(drop1, drop2);        // [J] surft / pi
    const double surf_c_pi = equivalent_surfenergy(drop1, drop2); // [J] surfc / pi
    const double etot_pi = kinetic_energy(drop1, drop2) + surf_t - surf_c; // [J] total energy / pi

    if (etot_pi < energylim)
    {
      const double exp = exponential(etot, surf_c);
      const double radiiratio = sizeratio_factor(drop1.radius,
                                                 drop2.radius);
      const double coaleff = aconst * radiiratio * exp;
      
      const double eff = colleff(drop1, drop2) * coaleff;
      return eff;
    }
    else // coaleff = 0.0
    {
      return 0.0;
    }
  }
};

HydrodynamicProb<LongKernelEfficiency, SimmelTerminalVelocity>
LongCollCoalProb()
/* returns the probability of collision-coalescence
using Simmel et al. 2002's formulation of
Long's Hydrodynamic Kernel */
{
  return HydrodynamicProb(LongKernelEfficiency{1.0},
                          SimmelTerminalVelocity{});
}

template <VelocityFormula TerminalVelocity>
HydrodynamicProb<LowListKernelEfficiency, TerminalVelocity>
LowListCoalProb(TerminalVelocity terminalv)
/* returns the probability of collision-coalescence
using Long's Hydrodynamic Kernel combined with
the coalescence efficiency from Low and List 1982. */
{
  return HydrodynamicProb(LowListKernelEfficiency{terminalv},
                          terminalv);
}

template <SDPairProbability CollisionXProbability>
SdmProcess auto CollisionCoalescenceProcess(const int interval,
                                            const std::function<double(int)> int2time,
                                            const CollisionXProbability p)
{
  const double realtstep = int2time(interval);

  CollisionX<CollisionXProbability, Coalescence>
      coal(realtstep, p, Coalescence{});

  return ConstTstepProcess{interval, coal};
}

#endif // COALESCENCE_HPP