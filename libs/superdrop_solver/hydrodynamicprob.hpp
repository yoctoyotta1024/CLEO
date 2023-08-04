// Author: Clara Bayley
// File: hydrodynamicprob.hpp
/* Header file for calculation of
probability of collision-coalescence
or collision-breakup between two
droplets using the hydrodynamic
ie. gravitational kernel */

#ifndef HYDRODYNAMICPROB_HPP
#define HYDRODYNAMICPROB_HPP

#include <concepts>
#include <cmath>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./terminalvelocity.hpp"

namespace dlc = dimless_constants;

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

#endif // HYDRODYNAMICPROB_HPP