// Author: Clara Bayley
// File: terminalvelocity.hpp
/* Header file for collision-coalescence 
probabilities for the SDM collision method
using various coalescence kernels.
Probs. are contained in structures to satisfy
requirements for using PairProbability 
concept in Collisionsmethod struct  */

#ifndef COALESCENCEKERNEL_HPP 
#define COALESCENCEKERNEL_HPP 

#include <algorithm>
#include <math.h>
#include <cmath>

#include "claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./terminalvelocity.hpp"

namespace dlc = dimless_constants

struct GolovinProb
{
  const double prob_jk_const;

  GolovinProb(const double R0)
      : prob_jk_const(1.5e3 * (pow(R0, 3.0))) {}

  const double operator()(const Superdrop &drop1,
                          const Superdrop &drop2,
                          const double DELT,
                          const double VOLUME) const
  /* returns probability that a pair of droplets coalesces
  according to Golovin's (sum of volumes) coalescence kernel. 
  Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol where 
  K(drop1, drop2) := C(drop1, drop2) * |v1−v2|, (see Shima 2009 eqn 3),
  and K(drop1, drop2) is Golovin's coalescence kernel */
  {
    const double DELT_DELVOL = DELT / VOLUME;                                   // time interval / volume for which collision probability is calculated [s/m^3]
    const double golovins_kernel = prob_jk_const * (drop1.vol() + drop2.vol()); // Golovin 1963 coalescence kernel
   
    const double prob_jk = golovins_kernel * DELT_DELVOL;
   
    return prob_jk;
  }
};

struct LongHydrodynamicProb
{
  const double prob_jk_const;
  const SimmelTerminalVelocity simmel_terminalv;

  LongHydrodynamicProb()
      : prob_jk_const(M_PI * pow(dlc::R0, 2.0) * dlc::W0),
        simmel_terminalv(SimmelTerminalVelocity{}) {}

  const double operator()(const Superdrop &drop1,
                          const Superdrop &drop2,
                          const double DELT,
                          const double VOLUME) const
  /* returns probability that a pair of droplets coalescese
  according to Long's (hydrodynamic i.e. gravitational) coalescence kernel. 
  Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol where 
  K(drop1, drop2) := C(drop1, drop2) * |v1−v2|, (see Shima 2009 eqn 3), 
  and K(drop1, drop2) is Long's coalescence kernel. Kernel equations taken
  from Simmel at al. 2002. colleff = collision efficiency E(R,r) where R>r.
  E(R,r) = E_coll(R,r) * E_coal(R,r) = E_coll(R,r) since it's assumed that
  E_coal(R,r) = 1. For collisions where R > rlim, E(R,r) = E_coll(R,r) = 1. */
  {
    constexpr double rlim = 5e-5 / dlc::R0; // 50 micron limit to determine collision efficiency (colleff)
    constexpr double smallcolleff_max = 0.001; // maximum collision efficiency if larger droplet's radius < rlim
    constexpr double A1 = 4.5e4 * dlc::R0 * dlc::R0; // constants in collision efficiency calc if larger droplet's radius < rlim
    constexpr double A2 = 3e-4 / dlc::R0;

    const double DELT_DELVOL = DELT / VOLUME;                                   // time interval / volume for which collision probability is calculated [s/m^3]
    const double bigr = std::max(drop1.radius, drop2.radius);
    const double smallr = std::min(drop1.radius, drop2.radius); 

    /* calculate collision efficiency */
    double colleff(1.0); // 
    if (bigr < rlim)
    {
      const double smallcolleff = A1 * pow(bigr, 2.0) * (1 - A2/smallr);
      colleff = std::max(smallcolleff, smallcolleff_max);
    }

    /* calculate Long's hydrodynamic (i.e. gravitational) 
    collision kernel according to Simmel et al. 2002 */
    const double v1 = simmel_terminalv(drop1);
    const double v2 = simmel_terminalv(drop2);
    const double longs_kernel = prob_jk_const * colleff * pow((bigr + smallr), 2.0) * std::abs(v1 - v2);
    
    const double prob_jk = longs_kernel * DELT_DELVOL;
    
    return prob_jk;
  }
};

#endif // COALESCENCEKERNEL_HPP 