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

struct LongHydrodynamicCoalProb
{
  const double prob_jk_const;
  const double coaleff;
  const SimmelTerminalVelocity terminalv;

  LongHydrodynamicCoalProb()
      : prob_jk_const(M_PI * pow(dlc::R0, 2.0) * dlc::W0),
        coaleff(1.0),
        terminalv(SimmelTerminalVelocity{}) {}

  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const
  /* returns probability that a pair of droplets coalescese
  according to Long's (hydrodynamic i.e. gravitational)
  collision-coalescence kernel. 
  Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol where 
  K(drop1, drop2) := C(drop1, drop2) * |v1−v2|, (see Shima 2009 eqn 3), 
  and K(drop1, drop2) is Long's collision-coalescence kernel.
  Kernel equations taken from Simmel at al. 2002.
  colleff = collision efficiency E(R,r) where R>r.
  E(R,r) = E_coll(R,r) * E_coal(R,r). Here it's assumed that
  E_coal(R,r) = 1. This also means for collisions where R > rlim, 
  E(R,r) = E_coll(R,r) = 1. */
  {
    constexpr double rlim = 5e-5 / dlc::R0; // 50 micron limit to determine collision-coalescence efficiency (eff)
    constexpr double smallcolleff_max = 0.001; // maximum efficiency if larger droplet's radius < rlim
    constexpr double A1 = 4.5e4 * dlc::R0 * dlc::R0; // constants in efficiency calc if larger droplet's radius < rlim
    constexpr double A2 = 3e-4 / dlc::R0;

    const double DELT_DELVOL = DELT / VOLUME;                                   // time interval / volume for which collision probability is calculated [s/m^3]
    const double bigr = std::max(drop1.radius, drop2.radius);
    const double smallr = std::min(drop1.radius, drop2.radius); 

    /* calculate collision-coalescence efficiency, eff = colleff * coaleff */
    double colleff(1.0); // 
    if (bigr < rlim)
    {
      const double smallcolleff = A1 * pow(bigr, 2.0) * (1 - A2/smallr);
      colleff = std::max(smallcolleff, smallcolleff_max);
    }
    const double eff = colleff * coaleff;

    /* calculate Long's hydrodynamic (i.e. gravitational) 
    collision-coalescence kernel according to Simmel et al. 2002 */
    const double sumrsqrd = pow((bigr + smallr), 2.0);
    const double vdiff = std::abs(terminalv(drop1) - terminalv(drop2));
    const double longs_kernel = prob_jk_const * eff * sumrsqrd * vdiff;
    
    const double prob_jk = longs_kernel * DELT_DELVOL;
    
    return prob_jk;
  }
};

#endif // COLLISIONKERNELS_HPP 