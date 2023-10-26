/*
 * ----- CLEO -----
 * File: impliciteuler.hpp
 * Project: superdrops
 * Created Date: Thursday 26th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 26th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header for class implementing the implicit euler method
 * for radial growth/shrink of each superdroplet due to
 * condensation and diffusion of water vapour according
 * to equations from "An Introduction To Clouds From
 * The Microscale to Climate" by Lohmann, Luond and
 * Mahrt, 1st edition." and Shima et al. 2009 */

#ifndef IMPLICITEULER_HPP
#define IMPLICITEULER_HPP

#include <math.h> // for fmax()

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "../cleoconstants.hpp"

namespace dlc = dimless_constants;

struct ImplicitIteration
{
  unsigned int niters; // number of NR iterations to try before testing convergence
  double subdelt;      // timestep of implicit method (at each substep >= niter NR iterations occur)
  double rtol;         // relative tolerance for convergence of NR method
  double atol;         // abolute tolerance for convergence of NR method

  double s_ratio;
  double akoh;
  double bkoh;
  double ffactor;

  KOKKOS_INLINE_FUNCTION double
  initialguess(const double rprev) const;
  /* returns appropriate initial value (ie. a reasonable guess) for
  'ziter' to use as first iteration of newton raphson method in
  rootfinding algorithm for timestepping condensation/evaporation ODE.
  Criteria is as in SCALE-SDM for making initial guess for given droplet
  much greater than its (activation radius)^2 if the
  supersaturation > its activation supersaturation  */

  KOKKOS_INLINE_FUNCTION
  double initialguess_shima(const double rprev) const;
  /* returns appropriate initial value (ie. a reasonable guess) for
  'ziter' to use as first iteration of newton raphson method in
  rootfinding algorithm for timestepping condensation/evaporation ODE.
  Criteria for modifying guess from rprev^2 are adapted from SCALE-SDM.
  Second criteria is that initial guess >= 'r1sqrd', where r1 is the
  equilibrium radius of a given droplet when s_ratio=1  */

  KOKKOS_INLINE_FUNCTION
  double newtonraphson_niterations(const double rprev,
                                   double ziter) const;
  /* Timestep condensation ODE by delt given initial guess for ziter,
  (which is usually radius^squared from previous timestep). Uses newton
  raphson iterative method to find new value of radius that converges
  on the root of the polynomial g(ziter) within the tolerances of the
  ImpIter instance. Uniquesol method assumes that solution to g(ziter)=0
  is unique and therefore Newton Raphson root finding algorithm converges
  quickly. This means method can be used with comparitively large tolerances
  and timesteps, and the maximum number of iterations is small. After
  'niters' iterations, convergence criteria is tested and futher
  iterations undertaken if not yet converged. */

  KOKKOS_INLINE_FUNCTION
  double newtonraphson_untilconverged(const unsigned int iterlimit,
                                      const double rprev,
                                      double ziter) const;
  /*  Timestep condensation ODE by delt given initial guess for ziter,
  (which is usually radius^squared from previous timestep). Uses
  newton raphson iterative method to find new value of radius that
  converges on the root of the polynomial g(ziter) within the tolerances
  of the ImpIter instance. After every iteration, convergence criteria is
  tested and error is raised if method does not converge within
  'iterlimit' iterations. Otherwise returns new value for the radius
  (which is the radius at timestep 't+subdelt'. Refer to section 5.1.2 Shima
  et al. 2009 and section 3.3.3 of Matsushima et al. 2023 for more details. */

  KOKKOS_INLINE_FUNCTION Kokkos::pair<bool, double>
  iterate_rootfinding_algorithm(const double rprev, double ziter) const;
  /* function performs one iteration of Newton Raphson rootfinding
  method and returns updated value of radius^2 alongside a boolean that
  is false if algorithm has converged */

  KOKKOS_INLINE_FUNCTION
  double ode_gfunc(const double rprev, const double rsqrd) const;
  /* returns g(z) / z * delt for g(z) function used in root finding
  Newton Raphson Method for dr/dt condensation / evaporation ODE.
  ODE is for radial growth/shrink of each superdroplet due to
  condensation and diffusion of water vapour according to
  equations from "An Introduction To Clouds...."
  (see note at top of file). Note: z = ziter = radius^2 */

  KOKKOS_INLINE_FUNCTION
  double ode_gfuncderivative(const double rsqrd) const;
  /* dg(z)/dz * delt, where dg(z)/dz is derivative of g(z) with
  respect to z=rsqrd. g(z) is polynomial to find root of using
  Newton Raphson Method. */

  KOKKOS_INLINE_FUNCTION
  bool isnotconverged(const double gfunciter,
                      const double gfuncprev) const
  /* boolean where True means criteria for ending newton raphson
  iterations has not yet been met. Criteria is standard local error
  test: |iteration - previous iteration| < RTOL * |iteration| + ATOL */
  {
    const double converged = rtol * Kokkos::abs(gfunciter) + atol;
    const double currentvalue = Kokkos::abs(gfunciter - gfuncprev);

    return (currentvalue >= converged); // true means it's not yet converged
  }
};

class ImplicitEuler
/* class for the Implicit Euler (IR) integration of
  to superdroplet condensational growth ODE. Uses a
  Newton Raphson (NR) root finding method to solve
  the implicit timestepping equation of stiff ODE */
{
private:
  unsigned int niters; // suggested number of iterations for implicit method before testing for convergence
  double delt;         // timestep of ODE solver (at each step implicit method is called)
  double maxrtol;      // adjustable relative tolerance for convergence of NR method
  double maxatol;      // adjustable abolute tolerance for convergence of NR method
  double subdelt;      // number of substeps to condensation method when supersat close to 1 (0 = false ie. none)

  KOKKOS_INLINE_FUNCTION
  double substepped_implicitmethod(const ImplicitIteration &implit,
                                   const unsigned int nsubsteps,
                                   const double rprev) const;

public:
  ImplicitEuler(const unsigned int niters,
                const double delt,
                const double maxrtol,
                const double maxatol,
                const double subdelt)
      : niters(niters), delt(delt), maxrtol(maxrtol),
        maxatol(maxatol), subdelt(subdelt) {}

  KOKKOS_INLINE_FUNCTION
  double solve_condensation_matsushima(const double s_ratio,
                                       const Kokkos::pair<double, double> akoh_bkoh,
                                       const double ffactor,
                                       const double rprev) const;
  /* forward timestep previous radius 'rprev' by delt using an implicit
  euler method to integrate the condensation/evaporation ODg. Implict
  timestepping equation defined in section 5.1.2 of Shima et al. 2009
  and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
  Newton Raphson iterations are used to converge towards the root of
  g(z) within the tolerances of an ImpIter instance. Tolerances,
  maxium number of iterations and sub-timestepping are adjusted based on
  on the uniqueness criteria of the polynomial g(z). Refer to section
  5.1.2 Shima et al. 2009 and section 3.3.3 of Matsushima et al. 2023
  for more details. */

  KOKKOS_INLINE_FUNCTION
  double solve_condensation(const double s_ratio,
                            const Kokkos::pair<double, double> akoh_bkoh,
                            const double ffactor,
                            const double rprev) const;
  /* forward timestep previous radius 'rprev' by delt using an implicit
  euler method to integrate the condensation/evaporation ODg. Implict
  timestepping equation defined in section 5.1.2 of Shima et al. 2009
  and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
  Newton Raphson iterations are used to converge towards the root of
  g(z) within the tolerances of an ImpIter instance. Tolerances,
  maxium number of iterations and sub-timestepping are adjusted when
  near to supersaturation=1 (when activation / deactivation may occur).
  Refer to section 5.1.2 Shima et al. 2009 and section 3.3.3 of
  Matsushima et al. 2023 for more details. */
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_INLINE_FUNCTION double
ImplicitEuler::solve_condensation(const double s_ratio,
                                  const Kokkos::pair<double, double> kohler_ab,
                                  const double ffactor,
                                  const double rprev) const
/* forward timestep previous radius 'rprev' by delt using an implicit
euler method to integrate the condensation/evaporation ODg. Implict
timestepping equation defined in section 5.1.2 of Shima et al. 2009
and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
Newton Raphson iterations are used to converge towards the root of
g(z) within the tolerances of an ImpIter instance. Tolerances,
maxium number of iterations and sub-timestepping are adjusted when
near to supersaturation=1 (when activation / deactivation may occur).
Refer to section 5.1.2 Shima et al. 2009 and section 3.3.3 of
Matsushima et al. 2023 for more details. */
{
  const double akoh(kohler_ab.first);
  const double bkoh(kohler_ab.second);

  const double s_act(1 + Kokkos::sqrt(4.0 * Kokkos::pow(akoh, 3.0) / 27 / bkoh)); // activation supersaturation
  if ((s_ratio > 0.999 * s_act) && (s_ratio < 1.001 * s_act))
  /* if supersaturation close to s_act, activation or
  deactivation might occur so perform subtimestepping */
  {
    const unsigned int nsubs = Kokkos::ceil(delt / subdelt);
    const double subt = delt / (double)nsubs;
    const ImplicitIteration implit{niters, subt, maxrtol, maxatol,
                                   s_ratio, akoh, bkoh, ffactor};

    return substepped_implicitmethod(implit, nsubs, rprev);
  }

  // else
  // /* Far from activation / deactivation appropriate choice
  // of initial guess allows rapid convergence of to correct
  // solution even in cases when spurious solutions exist
  // (see Matsushima et al. 2023 unqiuenss criteria). */
  // {
  //   const ImpIter impit{niters, delt, maxrtol, maxatol,
  //                       s_ratio, akoh, bkoh, ffactor};
  //   double init_ziter(impit.initialguess(rprev));
  //   return impit.newtonraphson_niterations(rprev, init_ziter);
  // }

  return rprev; // TODO delete
}

KOKKOS_INLINE_FUNCTION double
ImplicitEuler::substepped_implicitmethod(const ImplicitIteration &implit,
                                         const unsigned int nsubsteps,
                                         const double rprev) const
{
  double subr(rprev);
  for (unsigned int n = 0; n < nsubsteps; ++n)
  {
    double init_ziter(implit.initialguess(subr));
    // subr = implit.newtonraphson_niterations(subr, init_ziter);
  }
  return subr;
}

KOKKOS_INLINE_FUNCTION double
ImplicitIteration::initialguess(const double rprev) const
/* returns appropriate initial value (ie. a reasonable guess) for
'ziter' to use as first iteration of newton raphson method in
rootfinding algorithm for timestepping condensation/evaporation ODE.
Criteria is as in SCALE-SDM for making initial guess for given droplet
much greater than its (activation radius)^2 if the 
supersaturation > its activation supersaturation  */
{
  const double rprevsqrd(rprev * rprev);
  const double s_act(1 + Kokkos::sqrt(4.0 * Kokkos::pow(akoh, 3.0) / 27 / bkoh)); // activation supersaturation

  if (s_ratio > s_act)
  {
    constexpr double bigr(1e-3 / dlc::R0); // large initial guess for radius = 1mm for drop that should already be activated
    return fmax(bigr * bigr, rprevsqrd);   // Kokkos compatible equivalent to std::max() for floating point numbers
  }

  return rprevsqrd;
}

KOKKOS_INLINE_FUNCTION double
ImplicitIteration::initialguess_shima(const double rprev) const
/* returns appropriate initial value (ie. a reasonable guess) for
'ziter' to use as first iteration of newton raphson method in
rootfinding algorithm for timestepping condensation/evaporation ODE.
Criteria for modifying guess from rprev^2 are adapted from SCALE-SDM.
Second criteria is that initial guess >= 'r1sqrd', where r1 is the
equilibrium radius of a given droplet when s_ratio=1  */
{
  const double rsqrd(initialguess(rprev));
  const double r1sqrd(bkoh / akoh);
  return fmax(rsqrd, r1sqrd); // Kokkos compatible equivalent to std::max() for floating point numbers
}

#endif // IMPLICITEULER_HPP