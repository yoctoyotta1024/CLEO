/*
 * ----- CLEO -----
 * File: hydrodynamicprob.hpp
 * Project: collisionprobs
 * Created Date: Thursday 9th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 22nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for probability of some kind of
 * collision  event between two (real) droplets
 * using the hydrodynamic (i.e. gravitational)
 * kernel
*/

#ifndef HYDRODYNAMICPROB_HPP 
#define HYDRODYNAMICPROB_HPP 

#include <Kokkos_Core.hpp>

#include "../../cleoconstants.hpp"
#include "../superdrop.hpp"
#include "../terminalvelocity.hpp"


namespace dlc = dimless_constants;

template <typename E>
concept KernelEfficiency = requires(E e,
                                    const Superdrop &d1,
                                    const Superdrop &d2)
/* Objects that are of type 'KernelEfficiency'
take a pair of superdroplets and returns
something convertible to a double (such as the
efficiency factor for the hydrodynamic kernel) */
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
      : prob_jk_const(Kokkos::numbers::pi * dlc::R0 * dlc::R0 * dlc::W0),
        eff(e),
        terminalv(tv) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const
  /* returns probability that a pair of droplets collide
  (and coalesce or breakup etc.) according to the hydrodynamic,
  i.e. gravitational, collision kernel. Probability is given by
  prob_jk = K(drop1, drop2) * delta_t/delta_vol, (see Shima 2009 eqn 3)
  where the kernel, K(drop1, drop2) := eff * pi * (r1 + r2)^2 * |v1−v2|, 
  given the efficiency factor eff = eff(drop1, drop2), for
  example as expressed in equation 11 of Simmel at al. 2002 for
  collision-coalescence */
  {
    /* time interval / volume for which
    probability is calculated [s/m^3] */
    const double DELT_DELVOL = DELT / VOLUME;

    /* calculate Hydrodynamic Kernel*/
    const double sumr(drop1.get_radius() + drop2.get_radius());
    const double sumrsqrd(sumr * sumr);
    const double vdiff(Kokkos::abs(terminalv(drop1) - terminalv(drop2)));
    const double hydro_kernel(eff(drop1, drop2) *
                              prob_jk_const *
                              sumrsqrd * vdiff);

    /* calculate probability prob_jk analogous Shima 2009 eqn 3 */
    const double prob_jk = hydro_kernel * DELT_DELVOL;

    return prob_jk;
  }
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */


#endif // HYDRODYNAMICPROB_HPP 
