// Author: Clara Bayley
// File: collisionkernels.hpp
/* Header file for probabilities
of collision-coalescence or collision-breakup
for the SDM collisionx method using various
collision kernels. Probability calculations are
contained in structures that satisfy requirements
for using SDPairProbability
concept in CollisionX struct (see collisionx.hpp)  */

#ifndef COLLISIONKERNELS_HPP 
#define COLLISIONKERNELS_HPP 

#include <algorithm>
#include <math.h>
#include <cmath>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./terminalvelocity.hpp"

namespace dlc = dimless_constants;

struct GolovinCoalProb
{
  const double prob_jk_const;

  GolovinCoalProb(const double R0)
      : prob_jk_const(1.5e3 * (pow(R0, 3.0))) {}

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

template <typename E>
concept KernelEfficiency = requires(E e,
                                    const Superdrop &d1,
                                    const Superdrop &d2)
/* Objects that are of type 'KernelEfficiency'
take a pair of superdroplets and returns
something convertible to a double (such as the
efficiency factor for a collision kernel) */
{
  {
    e(d1, d2)
  } -> std::convertible_to<double>;
};

template <KernelEfficiency Efficiency,
          VelocityFormula TerminalVelocity>
struct HydrodynamicProb
{
  const double prob_jk_const;
  const Efficiency eff;
  const TerminalVelocity terminalv;

  HydrodynamicProb(Efficiency e, TerminalVelocity tv)
      : prob_jk_const(M_PI * pow(dlc::R0, 2.0) * dlc::W0),
        eff(e),
        terminalv(tv) {}

  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const
  /* returns probability that a pair of droplets collide (and coalesce
  or breakup) according to Long's formulation of the 
  hydrodynamic i.e. gravitational collision-interaction kernel.
  Probability equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol
  where K(drop1, drop2) := C(drop1, drop2) * |v1−v2|,
  (see Shima 2009 eqn 3), is the hydrodynamic collision-interaction
  kernel, for example expressed in equation 11 of Simmel at al. 2002
  for collision-coalescence */
  {
    /* time interval / volume for which
    probability is calculated [s/m^3] */
    const double DELT_DELVOL = DELT / VOLUME;  

    /* calculate Hydrodynamic Kernel*/                                 
    const double sumrsqrd = pow((drop1.radius + drop2.radius), 2.0);
    const double vdiff = std::abs(terminalv(drop1) - terminalv(drop2));
    const double hydro_kernel = prob_jk_const * sumrsqrd *
                                eff(drop1, drop2) * vdiff;

    /* calculate probability prob_jk analogous Shima 2009 eqn 3 */
    const double prob_jk = hydro_kernel * DELT_DELVOL;
    
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
                    const Superdrop &drop2)
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
      const double colleff_calc = A1 * pow(bigr, 2.0) * (1 - A2/smallr);
      colleff = std::max(colleff_calc, colleff_lim); // colleff >= colleff_lim
    }
    
    const double eff = colleff * coaleff;

    return eff;
  }
};

HydrodynamicProb<LongKernelEfficiency, SimmelTerminalVelocity>
LongHydrodynamicCollCoalProb()
{
  return HydrodynamicProb(LongKernelEfficiency{},
                          SimmelTerminalVelocity{});
}

#endif // COLLISIONKERNELS_HPP 