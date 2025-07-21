/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: impliciteuler.cpp
 * Project: superdrops
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors: Florian Poydenot
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

#include "impliciteuler.hpp"

/**
 * @brief Integrates the condensation / evaporation ODE employing the Implicit Euler method
 * similarly to Matsushima et. al, 2023.
 *
 * Forward timestep previous radius 'rprev' by delt using an Implicit Euler method (possibly with
 * sub-timestepping) to integrate the condensation/evaporation ODE using fixed thermodynamics from
 * the start of the timestep. Sub-timestepping employed when unique solution to g(Z) within
 * required radius range is not guarenteed. Criteria as in appendix C of Matsushima et. al, 2023
 * except minimum sub-timestep is limited by minsubdelt.
 *
 */
KOKKOS_FUNCTION double ImplicitEuler::solve_condensation(
    const double s_ratio, const Kokkos::pair<double, double> kohler_ab, const double ffactor,
    const double rprev) const {
  const auto ffactor_fv = ffactor / ventilation_factor(rprev);
  const auto odeconsts =
      ImplicitIterations::ODEConstants{s_ratio, kohler_ab.first, kohler_ab.second, ffactor_fv};

  auto ziter = implit.initialguess(odeconsts, rprev);
  const bool ucrit1 = first_unique_criteria(odeconsts, rprev, ziter);
  const bool ucrit2 = second_unique_criteria(odeconsts, delt);

  auto rsqrd = double{0.0};
  if (ucrit1 || ucrit2) {
    rsqrd = implit.integrate_condensation_ode(odeconsts, delt, rprev, ziter);
  } else {
    rsqrd = solve_with_adaptive_subtimestepping(odeconsts, delt, rprev, ziter);
  }

  return Kokkos::sqrt(rsqrd);
}

/**
 * @brief Test of uniqueness criteria for un-activated droplets in environment with
 * supersaturation less than its activation supersaturation.
 *
 * Returns true if solution to g(Z) is guarenteed to be unique because it meets the
 * uniquenes criteria of Case 2 from Matsushima et al. 2023 (see appendix C), namely that there
 * is only one real root to g(Z) in the range 0 < Z < critical_R^2, where critical_R is the
 * critical i.e. activation radius of the droplet. Here we use the less stringent constrain
 * that S <= S_crit rather than S <= 1, and we ensure the current value for ziter is also less
 * than the critical_R as it must be to guarentee solution in range 0 < R < critical_R
 * is converged upon.
 *
 */
KOKKOS_FUNCTION bool ImplicitEuler::first_unique_criteria(
    const ImplicitIterations::ODEConstants &odeconsts, const double rprev,
    const double ziter) const {
  const double akoh = odeconsts.akoh;
  const double bkoh = odeconsts.bkoh;

  const double rcritsqrd = 3.0 * bkoh / akoh;
  const bool is_ziter_unactivated = ziter < rcritsqrd;
  const bool is_unactivated = (rprev * rprev < rcritsqrd && is_ziter_unactivated);

  const double sqrdval = 4.0 * akoh * akoh * akoh / 27.0 / bkoh;
  const bool is_subactivated_saturation = (odeconsts.s_ratio <= 1.0 + Kokkos::pow(sqrdval, 0.5));

  return (is_unactivated && is_subactivated_saturation);
}

/**
 * @brief Integrates the condensation / evaporation ODE employing the Implicit Euler method
 * similarly to Matsushima et. al, 2023 with an adaptive timestepping subroutine.
 *
 * Forward timestep previous radius 'rprev' by delt using an Implicit Euler method with
 * sub-timestepping to integrate the condensation/evaporation ODE using fixed thermodynamics from
 * the start of the timestep. Sub-timestepping employed to try to ensure uique solution to g(Z)
 * as Matsushima et. al, 2023 except minimum sub-timestep is limited by minsubdelt.
 * If critdelt < minsubdelt the uniqueness is not guarenteed. Reducing minsubdelt therefore
 * increases the likelyhood of having a unique solution to g(Z), i.e. the accuracy of the
 * solver is increased.
 *
 */
KOKKOS_FUNCTION double ImplicitEuler::solve_with_adaptive_subtimestepping(
    const ImplicitIterations::ODEConstants &odeconsts, const double delt, double rprev,
    double ziter) const {
  const auto critdelt = critial_timestep(odeconsts);
  const auto mindelt = Kokkos::fmax(critdelt, minsubdelt);

  auto remdelt = delt;  // remaining time required to integrate over.
  while (remdelt > 0.0) {
    const auto subdelt = Kokkos::fmin(mindelt, remdelt);
    ziter = implit.integrate_condensation_ode(odeconsts, subdelt, rprev, ziter);
    rprev = Kokkos::pow(ziter, 0.5);
    remdelt -= subdelt;
  }

  return ziter;
}

/**
 * @brief Integrates the condensation / evaporation ODE for radius^2 from t -> t+ subdelt.
 *
 * Employs the Implicit Euler method (with potential sub-timestepping based on uniqueness criteria
 * of Matsushima et. al), 2023 to forward timestep previous radius 'rprev' by subdelt according to
 * the condensation/evaporation ODE. Implict timestepping equation defined in section 5.1.2 of
 * Shima et al. 2009 and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
 *
 * Uses at least 'niters' iterations of the Newton Raphson method and then checks if convergence
 * criteria has been met (if a root of the g(Z) polynomial has been converged upon), else performs
 * upto maxniters number of further iterations, checking for convergence after each one.
 *
 */
KOKKOS_FUNCTION
double ImplicitIterations::integrate_condensation_ode(const ODEConstants &odeconsts,
                                                      const double subdelt, const double rprev,
                                                      double ziter) const {
  const size_t niters = 2;
  const auto result =
      newtonraphson_niterations(odeconsts, subdelt, rprev, ziter, niters);  // ziter, is_converged

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
 */
KOKKOS_FUNCTION double ImplicitIterations::ode_gfunc(const ODEConstants &odeconsts,
                                                     const double subdelt, const double rprev,
                                                     const double rsqrd) const {
  const auto radius = double{Kokkos::sqrt(rsqrd)};

  const auto alpha = double{odeconsts.s_ratio - 1 - odeconsts.akoh / radius +
                            odeconsts.bkoh / Kokkos::pow(radius, 3.0)};
  const auto beta = double{2.0 * subdelt / (rsqrd * odeconsts.ffactor_fv)};
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
 */
KOKKOS_FUNCTION double ImplicitIterations::ode_gfuncderivative(const ODEConstants &odeconsts,
                                                               const double subdelt,
                                                               const double rsqrd) const {
  const auto radius = double{Kokkos::sqrt(rsqrd)};

  const auto alpha =
      double{odeconsts.akoh / radius - 3.0 * odeconsts.bkoh / Kokkos::pow(radius, 3.0)};
  const auto beta = double{subdelt / (rsqrd * odeconsts.ffactor_fv)};

  return 1 - alpha * beta;
}
