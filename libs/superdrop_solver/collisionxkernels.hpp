// Author: Clara Bayley
// File: collisionxkernels.hpp
/* Header file for calculation of probability
of collision-x (e.g collision-coalescence
or collusion-breakup) between two droplets
using a specific kernel e.g. Golovin's
or Long's or Low and List's. Probability
calculations are contained in structures
that satisfy the requirements of the
SDPairProbability concept also used by
the CollisionX struct */

#ifndef COLLISIONXKERNELS_HPP
#define COLLISIONXKERNELS_HPP

#include <concepts>
#include <cmath>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./terminalvelocity.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

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
private:
  const double prob_jk_const;
  const Efficiency eff;
  const TerminalVelocity terminalv;

public:
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

struct CollCoalProb_Golovin
/* Probability of collision-coalescence of
a pair of droplets according to Golovin 1963
(see e.g. Shima et al. 2009) */
{
private:
  const double prob_jk_const;

public:
  CollCoalProb_Golovin()
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

struct LongKernelEff
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

HydrodynamicProb<LongKernelEff, SimmelTerminalVelocity>
CollCoalProb_Long()
/* returns the probability of collision-coalescence
using Simmel et al. 2002's formulation of
Long's Hydrodynamic Kernel */
{
  return HydrodynamicProb(LongKernelEff{1.0},
                          SimmelTerminalVelocity{});
}


template <VelocityFormula TerminalVelocity>
struct CollisionKinetics
/* calculations involved in the kinetics of 
a collision between two superdroplets */
{
private:
  TerminalVelocity terminalv;

  /* constant required to calculate surface tension energy from
  dimensionless radius using surface tension of water = sigma = 7.28e-2 */
  const double surfconst{4.0 * 7.28e-2 * std::numbers::pi * dlc::R0 *dlc::R0}; // [J/m^-2]

public:
  CollisionKinetics(TerminalVelocity tv) : terminalv(tv) {};
  
  double collision_kinetic_energy(const Superdrop &drop1,
                                  const Superdrop &drop2) const
  /* returns cke/ pi, where cke = collision kinetic energy
  as formulated in Low and List 1982(a) eqn 3.1 */
  {
    constexpr double R0cubed = dlc::R0 * dlc::R0 * dlc::R0; // convert r^3 to [m^3]
    constexpr double ckeconst = R0cubed * 2.0 / 3.0 * DC::RHO_L *
                                std::numbers::pi * dlc::W0;
    
    const double r1_r2(drop1.radius / drop2.radius);
    const double rratio = std::pow(drop1.radius, 3.0) /
                          (1 + std::pow(r1_r2, 3.0)); // * R0cubed to convert to [m^3]

    const double vdiff = terminalv(drop1) - terminalv(drop2); // * dlc::W0 to convert to [m/s]

    const double cke = ckeconst * rratio * vdiff * vdiff;

    return cke;
  }

  double surfenergy(const Superdrop &drop) const
  /* returns energy due to surface tension of a single
  drop, analogous to equation 4.2 of Low and List 1982 */
  {
    const double rsqrd = drop.radius * drop.radius; // * R0sqrd to convert to [m^2]
    const double tot_surfe = surfconst * rsqrd;

    return tot_surfe; // total surface energy
  }

  double total_surfenergy(const Superdrop &drop1,
                          const Superdrop &drop2) const
  /* returns total energy due to surface tension of pair 
  of drops as in equation 4.2 of Low and List 1982 */
  {
    const double r1(drop1.radius);  
    const double r2(drop2.radius);  
    const double r2sum = (r1 * r1 + r2 * r2); // * R0sqrd to convert to [m^2]
    
    const double tot_surfe = surfconst * r2sum;

    return tot_surfe; // total surface energy
  }

  double coal_surfenergy(const Superdrop &drop1,
                         const Superdrop &drop2) const
  /* returns surface energy of single spherical equivalent, ie.
  coalesced state of two drops, divided by pi as in
  equation 4.3 of Low and List 1982 */
  {
    const double r1(drop1.radius); 
    const double r2(drop2.radius); 
    const double r3sum = std::pow(r1, 3.0) + std::pow(r2, 3.0);
    
    const double equiv_surfe = surfconst * std::pow(r3sum, 2.0 / 3.0);

    return equiv_surfe; // coalesced (spherical equivalent) surface energy
  }
};

template <VelocityFormula TerminalVelocity>
struct LowListCollCoalEff
/* coalescence and collision-coalescence efficiency factors
for the Hydrodynamic kernel. eff = colleff(R,r) * coaleff(R,r)
where:
- colleff is Long's collision efficiency as seen in
  equation 13 of Simmel et al. 2002
- coaleff is from equation (4.5) and (4.6) Low and List 1982(a) */
{
private:
  CollisionKinetics<TerminalVelocity> ck;
  LongKernelEff colleff{1.0};

  double exponential(const double etot,
                     const double surf_c) const
  /* calc exponential in eqn 4.5 Low and List 1982(a) given
  total collision energy, etot [J] and equivalent surface
  energy, surf_c [J] */
  {
    constexpr double bconst = -2.62e6; // [J^-2]
    constexpr double sigma = 7.28e-2;  // [J/m^-2]

    const double exponent = bconst * sigma * etot * etot / surf_c;

    return std::exp(exponent);
  }

  double sizeratio_factor(const double r1,
                          const double r2) const
  /* calc factor that takes into accoutn size ratio fo droplets in
  eqn 4.5 Low and List 1982(a). */
  {
    const double rsmall(std::min(r1, r2));
    const double rbig(std::max(r1, r2));
    const double alpha(1 + rsmall / rbig); // alpha = 1 + Ds/Dl

    return 1.0 / (alpha * alpha); // alpha^(-2)
  }

public:
  LowListCollCoalEff(TerminalVelocity tv) : ck(tv) {}

  double get_colleff(const Superdrop &drop1,
                     const Superdrop &drop2) const
  {
    return colleff(drop1, drop2);
  }               

  double coaleff(const Superdrop &drop1,
                 const Superdrop &drop2) const
  /* returns coaleff, the coalescence efficiency
  of two droplets (given that they have collided)
  from equation (4.5) and (4.6) Low and List 1982(a).
  the total collision-coealescence efficiency, 
  eff = coaleff * colleff, and the breakup efficiency 
  bueff = 1 - coaleff */
  {
    constexpr double aconst = 0.778;
    constexpr double energylim = 5e-6; // etot limit / pi [J]

    const double surf_t = ck.total_surfenergy(drop1, drop2);      // [J] surft / pi
    const double surf_c = ck.coal_surfenergy(drop1, drop2); // [J] surfc / pi
    const double etot = surf_t - surf_c +
                        ck.kinetic_energy(drop1, drop2); // [J] total energy / pi

    if (etot < energylim)
    {
      const double exp = exponential(etot, surf_c);
      const double radiiratio = sizeratio_factor(drop1.radius,
                                                 drop2.radius);
      const double coaleff = aconst * radiiratio * exp;

      return coaleff;
    }
    else // coaleff = 0.0
    {
      return 0.0;
    }
  }

  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2) const
  /* collision-coalescence efficiency, eff, using
  eff = colleff * coaleff */
  {
    const double eff = coaleff(drop1, drop2) * colleff(drop1, drop2); 

    return eff;
  }
};

template <VelocityFormula TerminalVelocity>
HydrodynamicProb<LowListCollCoalEff<TerminalVelocity>, TerminalVelocity>
CollCoalProb_LowList(TerminalVelocity tv)
/* returns the probability of collision-coalescence
using Long's Hydrodynamic Kernel combined with
the coalescence efficiency from Low and List 1982. */
{
  return HydrodynamicProb(LowListCollCoalEff(tv), tv);
}

template <VelocityFormula TerminalVelocity>
struct LowListCollBuEff
/* Collision-Breakup Efficiency factor, eff, for the
Hydrodynamic kernel. eff = colleff(R,r) * coaleff(R,r) where:
- colleff is Long's collision efficiency as seen in
  equation 13 of Simmel et al. 2002
- coaleff is from equation (4.5) and (4.6) Low and List 1982(a) */
{
private:
  LowListCollCoalEff<TerminalVelocity> lle;
  
public:
  LowListCollBuEff(TerminalVelocity tv) : lle(tv){};

  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2) const
  /* collision-breakup efficiency, eff, using
  eff = colleff * bueff, and bueff = 1 - coealeff
  as in McFarquhar 2004 (see equation (28) therein) */
  {
    const double bueff(1.0 - lle.coaleff(drop1, drop2));
    const double eff(bueff * lle.get_colleff(drop1, drop2) );

    return eff / 100.0; 
  }
};

template <VelocityFormula TerminalVelocity>
HydrodynamicProb<LowListCollBuEff<TerminalVelocity>, TerminalVelocity>
CollBuProb_LowList(TerminalVelocity tv)
/* returns the probability of collision-breakup
using Long's Hydrodynamic Kernel combined with
the breakup efficient from McFarquhar 2004 obtained from
the coalescence efficiency given by Low and List 1982. */
{
  return HydrodynamicProb(LowListCollBuEff(tv), tv);
}

#endif // COLLISIONXKERNELS_HPP