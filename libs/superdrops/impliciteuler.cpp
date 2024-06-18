/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: impliciteuler.cpp
 * Project: superdrops
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 18th June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for class implementing the Implicit Euler method for radial growth/shrink of each
 * droplet due to condensation / evaporation and diffusion of water vapour according
 * to equations from "An Introduction To Clouds From The Microscale to Climate" by Lohmann,
 * Luond and Mahrt, 1st edition." and Shima et al. 2009
 */

#include "./impliciteuler.hpp"

/**
 * @brief Integrates the condensation / evaporation ODE employing the Implicit Euler method
 * similarly to Matsushima et. al, 2023.
 *
 * Forward timestep previous radius 'rprev' by delt using an Implicit Euler method (possibly with
 * sub-timestepping) to integrate the condensation/evaporation ODE using fixed thermodynamics from
 * the start of the timestep.
 *
 * Sub-timestepping is used to ensure solution to g(Z) is unique, unless timestep is less than
 * minsubdelt.  TODO(CB): WIP <-
 *
 * @param s_ratio The saturation ratio.
 * @param kohler_ab A pair containing 'a' and 'b' factors for Kohler curve in that order.
 * @param ffactor The sum of the diffusion factors.
 * @param rprev previous radius at time = t
 * @return updated radius for time = t + delt
 */
KOKKOS_FUNCTION double ImplicitEuler::solve_condensation(
    const double s_ratio, const Kokkos::pair<double, double> kohler_ab, const double ffactor,
    const double rprev) const {
  const auto odeconsts =
      ImplicitIterations::ODEConstants{s_ratio, kohler_ab.first, kohler_ab.second, ffactor};
  const auto subdelt = delt;
  const auto rsqrd = implit.integrate_condensation_ode(odeconsts, subdelt, rprev);
  return Kokkos::sqrt(rsqrd);
}

// /**
//  *  TODO(CB): WIP ->
//  *
//  * @brief Performs the implicit method with sub-stepping.
//  *
//  * This method performs the implicit method with substepping, iterating over the substeps to
//  * compute the implicit method for each substep.
//  *
//  * @param implit object defining implicit iterations.
//  * @param nsubsteps number of substeps to perform.
//  * @param rprev Radius of droplet at previous timestep.
//  * @return New value for droplet radius.
//  */
// KOKKOS_FUNCTION double substepped_implicitmethod(const ImplicitIteration &implit,
//                                                                 const unsigned int nsubsteps,
//                                                                 const double rprev) const {
//   auto subr = rprev;
//   for (unsigned int n(0); n < nsubsteps; ++n) {
//     auto init_ziter = double{implit.initialguess(subr)};
//     subr = implit.newtonraphson_niterations(subr, init_ziter);
//   }
//   return subr;
// }

/**
 * @brief Integrates the condensation / evaporation ODE for radius^2 from t -> t+ subdelt.
 *
 * Employs the Implicit Euler method (with potential sub-timestepping based on uniqueness criteria
 * of Matsushima et. al), 2023 to forward timestep previous radius 'rprev' by subdelt according to
 * the condensation/evaporation ODE. Implict timestepping equation defined in section 5.1.2 of
 * Shima et al. 2009 and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
 *
 * Uses at least 2 iterations of the Newton Raphson method and then checks if convergence
 * criteria has been met (if a root of the g(Z) polynomial has been converged upon), else performs
 * upto maxniters number of further iterations, checking for convergence after each one.
 *
 * @param odeconsts Constants of ODE during integration
 * @param subdelt time over which to integrate ODE
 * @param rprev Radius of droplet at previous timestep.
 */
KOKKOS_FUNCTION double ImplicitIterations::integrate_condensation_ode(const ODEConstants &odeconsts,
                                                                      const double subdelt,
                                                                      const double rprev) const {
  auto ziter = initialguess(odeconsts, rprev);
  const auto result =
      newtonraphson_niterations(odeconsts, subdelt, rprev, ziter, 2);  // ziter, is_converged

  if (result.second) {
    return result.first;
  } else {
    return newtonraphson_untilconverged(odeconsts, maxniters, subdelt, rprev, ziter);
  }
}

/**
 * @brief Returns appropriate initial guess (ie. a reasonable guess) for the Newton-Raphson
 * method.
 *
 * This method returns an initial guess for the Newton-Raphson method based on
 * the given radius from the previous timestep and the current supersaturation ratio.
 *
 * Guess is supposed to be a reasonable value for initial 'ziter' to use as first iteration of NR
 * method in rootfinding algorithm for timestepping condensation/evaporation ODE. Here the guess
 * criteria are as in SCALE-SDM for making initial guess for given droplet much greater
 * than its (activation radius)^2 if the supersaturation > its activation supersaturation.
 *
 * @param odeconsts Constants of ODE during integration
 * @param rprev Radius of droplet at previous timestep.
 * @return Initial guess for ziter.
 */
KOKKOS_FUNCTION double ImplicitIterations::initialguess(const ODEConstants &odeconsts,
                                                        const double rprev) const {
  const auto s_act = double{1 + Kokkos::sqrt(4.0 * Kokkos::pow(odeconsts.akoh, 3.0) / 27 /
                                             odeconsts.bkoh)};  // activation supersaturation

  if (odeconsts.s_ratio > s_act) {
    constexpr double bigr(
        1e-3 /
        dlc::R0);  // large initial guess for radius = 1mm for drop that should already be activated
    const auto maxradius = Kokkos::fmax(bigr, rprev);
    return maxradius * maxradius;
  }

  return rprev * rprev;
}

/**
 * @brief Performs niters number of Newton-Raphson iterations.
 *
 * Function integrates (timesteps) condensation ODE by delt given initial guess for ziter,
 * (which is usually radius^squared from previous timestep). Uses Newton Raphson iterative method
 * with 'niters' number of iterations then returns updated ziter and boolean which is true if
 * rootfinding has passed convergence test.
 *
 * @param odeconsts Constants of ODE during integration
 * @param subdelt Time over which to integrate ODE
 * @param rprev Radius at previous timestep
 * @param ziter The current guess for ziter.
 * @param niters Number of iterations of NR method to perform
 * @return The updated value of ziter.
 */
KOKKOS_FUNCTION Kokkos::pair<double, bool> ImplicitIterations::newtonraphson_niterations(
    const ODEConstants &odeconsts, const double subdelt, const double rprev, double ziter,
    const size_t niters) const {
  auto is_converged = false;

  for (size_t iter(0); iter < niters; ++iter) {
    const auto result =
        iterate_rootfinding_algorithm(odeconsts, subdelt, rprev, ziter);  // ziter, is_converged
    ziter = result.first;
    is_converged = result.second;
  }

  return {ziter, is_converged};
}

/**
 *
 * @brief Performs Newton-Raphson iterations until convergence or maximum number of
 * iterations is reached.
 *
 * After every iteration, convergence criteria is tested and error is raised if method does not
 * converge within 'maxniters' iterations. Otherwise once convergence test is passed, function
 * returns the new value for the ziter (which is the radius^2 at timestep 't+delt'). Refer to
 * section 5.1.2 Shima et al. 2009 and section 3.3.3 of Matsushima et al. 2023 for more details.
 *
 * @param odeconsts Constants of ODE during integration
 * @param subdelt Time over which to integrate ODE
 * @param niterslimit The maxiumum number of iterations to attempt.
 * @param rprev Radius at the previous timestep.
 * @param ziter The current guess for ziter.
 * @return The updated value of ziter.
 */
KOKKOS_FUNCTION double ImplicitIterations::newtonraphson_untilconverged(
    const ODEConstants &odeconsts, const size_t niterslimit, const double subdelt,
    const double rprev, double ziter) const {
  auto is_converged = bool{false};
  auto niter = size_t{1};

  while (!is_converged) {
    assert((niter <= niterslimit) &&
           "No root converged upon within max number of iterations of Newton Raphson Method.");
    const auto result =
        iterate_rootfinding_algorithm(odeconsts, subdelt, rprev, ziter);  // ziter, is_converged
    ziter = result.first;
    is_converged = result.second;
    niter += 1;
  }

  return ziter;
}

/**
 * @brief Perform one iteration of the Newton-Raphson rootfinding algorithm.
 *
 * This function performs one iteration of the Newton-Raphson rootfinding algorithm, i.e.
 * ziter^(m) -> ziter^(m+1) for iteration m+1 starting at m=1. Returns the updated value of ziter
 * alongside a boolean which is true if the new value of ziter passes the convergence test.
 *
 * _Note:_ ziter is limited to >= 1e-8 so it's always > 0.0
 *
 * @param odeconsts Constants of ODE during integration
 * @param subdelt Time over which to integrate ODE
 * @param rprev Radius at the previous timestep.
 * @param ziter The current guess for ziter.
 * @return A pair of the updated ziter and a boolean which is true if a root is converged upon.
 */
KOKKOS_FUNCTION
Kokkos::pair<double, bool> ImplicitIterations::iterate_rootfinding_algorithm(
    const ODEConstants &odeconsts, const double subdelt, const double rprev, double ziter) const {
  // perform iteration
  const auto numerator = ode_gfunc(odeconsts, subdelt, rprev, ziter);
  const auto denominator = ode_gfuncderivative(odeconsts, subdelt, ziter);
  ziter = ziter * (1 - numerator / denominator);

  // ensure ziter > 0.0
  ziter = Kokkos::fmax(ziter, 1e-8);

  // check if root has been converged upon
  const double newnumerator = ode_gfunc(odeconsts, subdelt, rprev, ziter);
  const auto is_converged = check_for_convergence(newnumerator, numerator);

  return {ziter, is_converged};
}

/**
 * @brief Returns the value of g(z) / z * subdelt for the ODE.
 *
 * This method computes the value of g(z) / z * subdelt used in the root-finding
 * Newton-Raphson method for the dr/dt condensation / evaporation ODE.
 *
 * ODE is for radial growth/shrink of each superdroplet due to condensation and diffusion
 * of water vapour according to equations from "An Introduction To Clouds...." (see note at
 * top of file).
 *
 * _Note:_ z = ziter = radius^2.
 *
 * @param odeconsts Constants of ODE during integration
 * @param subdelt Time over which to integrate ODE
 * @param rprev Radius at the previous timestep.
 * @param rsqrd Current radius squared.
 * @param subdelt Change in time to forward integrate ODE over.
 * @return RHS of g(z) / z * subdelt evaluted at rqrd.
 */
KOKKOS_FUNCTION double ImplicitIterations::ode_gfunc(const ODEConstants &odeconsts,
                                                     const double subdelt, const double rprev,
                                                     const double rsqrd) const {
  const auto radius = double{Kokkos::sqrt(rsqrd)};

  const auto alpha = double{odeconsts.s_ratio - 1 - odeconsts.akoh / radius +
                            odeconsts.bkoh / Kokkos::pow(radius, 3.0)};
  const auto beta = double{2.0 * subdelt / (rsqrd * odeconsts.ffactor)};
  const auto gamma = double{Kokkos::pow(rprev / radius, 2.0)};

  return 1 - gamma - alpha * beta;
}

/**
 * @brief Returns the value of the derivative of g(z) with respect to z.
 *
 * This method computes the value of dg(z)/dz * subdelt, where dg(z)/dz is the derivative of
 * g(z) with respect to z=rsqr. g(z) is polynomial to find root of using Newton Raphson Method
 * consistent as in ode_gfunc(...).
 *
 * @param odeconsts Constants of ODE during integration
 * @param subdelt Time over which to integrate ODE
 * @param rsqrd Current radius squared.
 * @param subdelt Change in time to forward integrate ODE over.
 * @return RHS of dg(z)/dz * subdelt evaluted at rqrd.
 */
KOKKOS_FUNCTION double ImplicitIterations::ode_gfuncderivative(const ODEConstants &odeconsts,
                                                               const double subdelt,
                                                               const double rsqrd) const {
  const auto radius = double{Kokkos::sqrt(rsqrd)};

  const auto alpha =
      double{odeconsts.akoh / radius - 3.0 * odeconsts.bkoh / Kokkos::pow(radius, 3.0)};
  const auto beta = double{subdelt / (rsqrd * odeconsts.ffactor)};

  return 1 - alpha * beta;
}

// /**
//  * TODO(CB): WIP ->
//  *
//  * @brief Integrates the condensation / evaporation ODE employing the Implicit Euler method
//  * as in Matsushima et al, 2023.
//  *
//  * Forward timestep previous radius 'rprev' by delt using an Implicit Euler method to integrate
//  * the condensation/evaporation ODE. Implict timestepping equation defined in section 5.1.2 of
//  * Shima et al. 2009 and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
//  *
//  * Newton Raphson iterations are used to converge towards the root of g(z) within the tolerances
//  * of an ImpIter instance. Tolerances, maxium number of iterations and sub-timestepping are
//  * adjusted based on the uniqueness criteria of the polynomial g(z). Uniqueness criteria, ucrit1
//  * and / or ucrit2, assume that solution to g(ziter)=0 is unique and therefore Newton Raphson
//  root
//  * finding algorithm converges quickly. This means method can be used with comparitively large
//  * tolerances and timesteps, and the maximum number of iterations is small. Refer to
//  section 5.1.2
//  * of Shima et al. 2009 and section 3.3.3 of Matsushima et al. 2023 for more details.
//  */
// KOKKOS_FUNCTION double ImplicitEuler::solve_condensation_matsushima(
//     const double s_ratio, const Kokkos::pair<double, double> kohler_ab, const double ffactor,
//     const double rprev) const {
//   const auto akoh = double{kohler_ab.first};
//   const auto bkoh = double{kohler_ab.second};

//   const auto max_uniquedelt = double{2.5 * ffactor / akoh * Kokkos::pow(5.0 * bkoh / akoh, 1.5)};
//   const auto ract_ratio = double{rprev * rprev * akoh / 3.0 / bkoh};

//   const auto ucrit1 = bool{(s_ratio <= 1.0 && ract_ratio < 1.0)};
//   const auto ucrit2 = bool{delt <= max_uniquedelt};

//   /* at least one criteria is met such that there is unique solution */
//   if (ucrit1 || ucrit2) {
//     const ImplicitIteration implit{niters, delt, maxrtol, maxatol, s_ratio, akoh, bkoh, ffactor};
//     auto init_ziter = double{implit.initialguess(rprev)};
//     return implit.newtonraphson_niterations(rprev, init_ziter);

//     /* In general there may be > 0 spurious solutions.
//     Convergence may be slower so allow >= 3 Newton Raphson
//     iterations (could also refine tolerances) */
//   } else {
//     auto subt = double{Kokkos::fmax(
//         max_uniquedelt,
//         subdelt)};  // Kokkos compatible equivalent to std::max() for floating point numbers
//     const auto nsubs = (unsigned int)Kokkos::ceil(delt / subt);
//     subt = delt / static_cast<double>(nsubs);

//     const ImplicitIteration implit{niters, subt, maxrtol, maxatol, s_ratio, akoh, bkoh, ffactor};

//     return substepped_implicitmethod(implit, nsubs, rprev);
//   }
// }
