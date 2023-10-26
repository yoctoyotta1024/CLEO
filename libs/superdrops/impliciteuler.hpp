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

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

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
  double substepped_implicitmethod(const unsigned int iters,
                                   const double delt,
                                   const double rtol,
                                   const double atol,
                                   const double s_ratio,
                                   const double akoh,
                                   const double bkoh,
                                   const double ffactor,
                                   const double rprev,
                                   double subdelt) const;

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

KOKKOS_INLINE_FUNCTION
double ImplicitEuler::solve_condensation(const double s_ratio,
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

  const double s_act(1 + std::sqrt(4.0 * std::pow(akoh, 3.0) / 27 / bkoh)); // activation supersaturation
  if ((s_ratio > 0.999 * s_act) && (s_ratio < 1.001 * s_act))
  /* if supersaturation close to s_act, activation or
  deactivation might occur so perform subtimestepping */
  {
    return substepped_implicitmethod(niters, delt, maxrtol, maxatol,
                                     s_ratio, akoh, bkoh, ffactor,
                                     rprev, subdelt);
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

KOKKOS_INLINE_FUNCTION
double ImplicitEuler::substepped_implicitmethod(const unsigned int iters,
                                                const double delt,
                                                const double rtol,
                                                const double atol,
                                                const double s_ratio,
                                                const double akoh,
                                                const double bkoh,
                                                const double ffactor,
                                                const double rprev,
                                                double subdelt) const
{
  const unsigned int nsubs = std::ceil(delt / subdelt);
  subdelt = delt / (double)nsubs;

  // const ImpIter impit{niters, subdelt, maxrtol, maxatol,
  //                       s_ratio, akoh, bkoh, ffactor};

  double subr(rprev);
  // for (unsigned int n=0; n < nsubs; ++n)
  // {
  //   double init_ziter(impit.initialguess(subr));
  //   subr = impit.newtonraphson_niterations(subr, init_ziter);
  // }
  return subr;
}

#endif // IMPLICITEULER_HPP