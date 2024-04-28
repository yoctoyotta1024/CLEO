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
 * Last Modified: Sunday 21st April 2024
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
 * as in Shima et. al, 2023 with adjustments for near-supersaturation conditions.
 *
 * Forward timestep previous radius 'rprev' by delt using an Implicit Euler method to integrate
 * the condensation/evaporation ODE. Implict timestepping equation defined in section 5.1.2 of
 * Shima et al. 2009 and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
 *
 * Newton Raphson iterations are used to converge towards the root of g(z) within the tolerances
 * of an ImpIter instance. Tolerances, maxium number of iterations and sub-timestepping are
 * adjusted when near to supersaturation=1 (when activation / deactivation may occur). Far from
 * activation, solution to g(ziter)=0 is usually unique and Newton Raphson root finding algorithm
 * converges quickly. This means method can be used with comparitively large tolerances and
 * timesteps, and the maximum number of iterations is small.
 */
KOKKOS_FUNCTION double ImplicitEuler::solve_condensation(
    const double s_ratio, const Kokkos::pair<double, double> kohler_ab, const double ffactor,
    const double rprev) const {
  const auto akoh = double{kohler_ab.first};
  const auto bkoh = double{kohler_ab.second};

  const auto s_act = double{
      1 + Kokkos::sqrt(4.0 / 27.0 * Kokkos::pow(akoh, 3.0) / bkoh)};  // activation supersaturation

  /* if supersaturation close to s_act, activation or
  deactivation might occur so perform subtimestepping */
  if ((s_ratio > 0.999 * s_act) && (s_ratio < 1.001 * s_act)) {
    const auto nsubs = (unsigned int)Kokkos::ceil(delt / subdelt);
    const auto subt = double{delt / static_cast<double>(nsubs)};
    const ImplicitIteration implit{niters, subt, maxrtol, maxatol, s_ratio, akoh, bkoh, ffactor};

    return substepped_implicitmethod(implit, nsubs, rprev);

    /* Far from activation / deactivation appropriate choice
    of initial guess allows rapid convergence of to correct
    solution even in cases when spurious solutions exist
    (see Matsushima et al. 2023 unqiuenss criteria). */
  } else {
    const ImplicitIteration implit{niters, delt, maxrtol, maxatol, s_ratio, akoh, bkoh, ffactor};
    auto init_ziter = double{implit.initialguess(rprev)};
    return implit.newtonraphson_niterations(rprev, init_ziter);
  }
}

/**
 * @brief Integrates the condensation / evaporation ODE employing the Implicit Euler method
 * as in Matsushima et al, 2023.
 *
 * Forward timestep previous radius 'rprev' by delt using an Implicit Euler method to integrate
 * the condensation/evaporation ODE. Implict timestepping equation defined in section 5.1.2 of
 * Shima et al. 2009 and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
 *
 * Newton Raphson iterations are used to converge towards the root of g(z) within the tolerances
 * of an ImpIter instance. Tolerances, maxium number of iterations and sub-timestepping are
 * adjusted based on the uniqueness criteria of the polynomial g(z). Uniqueness criteria, ucrit1
 * and / or ucrit2, assume that solution to g(ziter)=0 is unique and therefore Newton Raphson root
 * finding algorithm converges quickly. This means method can be used with comparitively large
 * tolerances and timesteps, and the maximum number of iterations is small. Refer to section 5.1.2
 * of Shima et al. 2009 and section 3.3.3 of Matsushima et al. 2023 for more details.
 */
KOKKOS_FUNCTION double ImplicitEuler::solve_condensation_matsushima(
    const double s_ratio, const Kokkos::pair<double, double> kohler_ab, const double ffactor,
    const double rprev) const {
  const auto akoh = double{kohler_ab.first};
  const auto bkoh = double{kohler_ab.second};

  const auto max_uniquedelt = double{2.5 * ffactor / akoh * Kokkos::pow(5.0 * bkoh / akoh, 1.5)};
  const auto ract_ratio = double{rprev * rprev * akoh / 3.0 / bkoh};

  const auto ucrit1 = bool{(s_ratio <= 1.0 && ract_ratio < 1.0)};
  const auto ucrit2 = bool{delt <= max_uniquedelt};

  /* at least one criteria is met such that there is unique solution */
  if (ucrit1 || ucrit2) {
    const ImplicitIteration implit{niters, delt, maxrtol, maxatol, s_ratio, akoh, bkoh, ffactor};
    auto init_ziter = double{implit.initialguess(rprev)};
    return implit.newtonraphson_niterations(rprev, init_ziter);

    /* In general there may be > 0 spurious solutions.
    Convergence may be slower so allow >= 3 Newton Raphson
    iterations (could also refine tolerances) */
  } else {
    auto subt = double{Kokkos::fmax(
        max_uniquedelt,
        subdelt)};  // Kokkos compatible equivalent to std::max() for floating point numbers
    const auto nsubs = (unsigned int)Kokkos::ceil(delt / subt);
    subt = delt / static_cast<double>(nsubs);

    const ImplicitIteration implit{niters, subt, maxrtol, maxatol, s_ratio, akoh, bkoh, ffactor};

    return substepped_implicitmethod(implit, nsubs, rprev);
  }
}

/**
 * @brief Performs the implicit method with sub-stepping.
 *
 * This method performs the implicit method with substepping, iterating over the substeps to
 * compute the implicit method for each substep.
 *
 * @param implit object defining implicit iterations.
 * @param nsubsteps number of substeps to perform.
 * @param rprev Radius of droplet at previous timestep.
 * @return New value for droplet radius.
 */
KOKKOS_FUNCTION double ImplicitEuler::substepped_implicitmethod(const ImplicitIteration &implit,
                                                                const unsigned int nsubsteps,
                                                                const double rprev) const {
  auto subr = rprev;
  for (unsigned int n(0); n < nsubsteps; ++n) {
    auto init_ziter = double{implit.initialguess(subr)};
    subr = implit.newtonraphson_niterations(subr, init_ziter);
  }
  return subr;
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
 * @param rprev Radius of droplet at previous timestep.
 * @return Initial guess for ziter.
 */
KOKKOS_FUNCTION double ImplicitIteration::initialguess(const double rprev) const {
  const auto rprevsqrd = double{rprev * rprev};
  const auto s_act = double{
      1 + Kokkos::sqrt(4.0 * Kokkos::pow(akoh, 3.0) / 27 / bkoh)};  // activation supersaturation

  if (s_ratio > s_act) {
    constexpr double bigr(
        1e-3 /
        dlc::R0);  // large initial guess for radius = 1mm for drop that should already be activated
    return Kokkos::fmax(
        bigr * bigr,
        rprevsqrd);  // Kokkos compatible equivalent to std::max() for floating point numbers
  }

  return rprevsqrd;
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
 * criteria are adapted from SCALE-SDM. Second criteria is that initial guess >= 'r1sqrd', where
 * r1 is the equilibrium radius of a given droplet when s_ratio=1.
 *
 * @param rprev Radius of droplet at previous timestep.
 * @return Initial guess for ziter.
 */
KOKKOS_FUNCTION double ImplicitIteration::initialguess_shima(const double rprev) const {
  const auto rsqrd = double{initialguess(rprev)};
  const auto r1sqrd = double{bkoh / akoh};
  return Kokkos::fmax(
      rsqrd, r1sqrd);  // Kokkos compatible equivalent to std::max() for floating point numbers
}

/**
 * @brief Performs Newton-Raphson iterations with a fixed number of iterations.
 *
 * This method performs Newton-Raphson iterations for a fixed number of iterations,
 * then tests for convergence and continues with futher iterations until convergence occurs or
 * or until the maximum number of iterations is reached.
 *
 * Funciton intergrates (timesteps) condensation ODE by delt given initial guess for ziter,
 * (which is usually radius^squared from previous timestep). Uses Newton Raphson iterative method
 * to find new value of the radius that converges on the root of the polynomial g(ziter) within
 * the tolerances of the ImpIter instance. After 'niters' iterations, convergence criteria is
 * tested and futher iterations undertaken if polynomial root has not yet been converged upon.
 *
 * @param rprev Radius at previous timestep
 * @param ziter The current guess for ziter.
 * @return The updated value of ziter.
 */
KOKKOS_FUNCTION double ImplicitIteration::newtonraphson_niterations(const double rprev,
                                                                    double ziter) const {
  // perform 'niters' iterations
  auto numerator = double{0.0};
  for (unsigned int iter(0); iter < niters; ++iter) {
    /* perform one attempted iteration  ziter^(m) -> ziter^(m+1)
    for iteration m+1 starting at m=1 */
    numerator = ode_gfunc(rprev, ziter);
    const auto denominator = ode_gfuncderivative(ziter);
    ziter -= ziter * numerator / denominator;  // increment ziter
    ziter = Kokkos::fmax(ziter, 1e-8);         // do not allow ziter < 0.0
  }

  // perform upto 'iterlimit' further iterations if convergence test fails
  if (!(isnot_converged(ode_gfunc(rprev, ziter), numerator))) {
    return Kokkos::sqrt(ziter);
  } else {
    constexpr unsigned int iterlimit = 50;  // maximum number of further iterations
    return newtonraphson_untilconverged(iterlimit, rprev, ziter);
  }
}

/**
 * @brief Returns the value of g(z) / z * delt for the ODE.
 *
 * This method computes the value of g(z) / z * delt used in the root-finding
 * Newton-Raphson method for the dr/dt condensation / evaporation ODE.
 *
 * ODE is for radial growth/shrink of each superdroplet due to condensation and diffusion
 * of water vapour according to equations from "An Introduction To Clouds...." (see note at
 * top of file).
 *
 * _Note:_ z = ziter = radius^2.
 *
 * @param rprev Radius at the previous timestep.
 * @param rsqrd Current radius squared.
 * @return RHS of g(z) / z * delt evaluted at rqrd.
 */
KOKKOS_FUNCTION double ImplicitIteration::ode_gfunc(const double rprev, const double rsqrd) const {
  const auto radius = double{Kokkos::sqrt(rsqrd)};

  const auto alpha = double{s_ratio - 1 - akoh / radius + bkoh / Kokkos::pow(radius, 3.0)};
  const auto beta = double{2.0 * subdelt / (rsqrd * ffactor)};
  const auto gamma = double{Kokkos::pow(rprev / radius, 2.0)};

  return 1 - gamma - alpha * beta;
}

/**
 * @brief Returns the value of the derivative of g(z) with respect to z.
 *
 * This method computes the value of dg(z)/dz * delt, where dg(z)/dz is the derivative of
 * g(z) with respect to z=rsqr. g(z) is polynomial to find root of using Newton Raphson Method
 * consistent as in ode_gfunc(...).
 *
 * @param rsqrd Current radius squared.
 * @return RHS of dg(z)/dz * delt evaluted at rqrd.
 */
KOKKOS_FUNCTION double ImplicitIteration::ode_gfuncderivative(const double rsqrd) const {
  const auto radius = double{Kokkos::sqrt(rsqrd)};

  const auto alpha = double{akoh / radius - 3.0 * bkoh / Kokkos::pow(radius, 3.0)};
  const auto beta = double{subdelt / (rsqrd * ffactor)};

  return 1 - alpha * beta;
}

/**
 * @brief Performs Newton-Raphson iterations until convergence.
 *
 * After every iteration, convergence criteria is tested and error is raised if method does not
 * converge within 'iterlimit' iterations. Otherwise returns new value for the radius (which is
 * the radius at timestep 't+subdelt'. Refer to section 5.1.2 Shima et al. 2009 and
 * section 3.3.3 of Matsushima et al. 2023 for more details.
 *
 * @param iterlimit The maxiumum number of iterations to attempt.
 * @param rprev Radius at the previous timestep.
 * @param ziter The current guess for ziter.
 * @return The updated value of ziter.
 */
KOKKOS_FUNCTION double ImplicitIteration::newtonraphson_untilconverged(const unsigned int iterlimit,
                                                                       const double rprev,
                                                                       double ziter) const {
  auto do_iter = bool{true};
  auto iter = (unsigned int)1;

  // perform newton raphson iterations if convergence test fails
  // and throw error if not converged within 'iterlimit' iterations
  while (do_iter) {
    assert((iter <= iterlimit) && "Newton Raphson Method not converged.");

    /* perform one attempted iteration  ziter^(m) -> ziter^(m+1)
    for iteration m+1 starting at m=1 and then test for convergence */
    const auto iterreturn(iterate_rootfinding_algorithm(rprev, ziter));
    do_iter = iterreturn.first;
    ziter = Kokkos::fmax(iterreturn.second, 1e-8);  // do not allow ziter < 0.0 (fmax ~ std::max())
    iter += 1;
  }

  return Kokkos::sqrt(ziter);
}

/**
 * @brief Perform one iteration of the Newton-Raphson rootfinding algorithm.
 *
 * This function performs one iteration of the Newton-Raphson rootfinding algorithm and returns
 * the updated value of radius^2 alongside a boolean indicating whether the algorithm has
 * converged.
 *
 * @param rprev Radius at the previous timestep.
 * @param ziter The current guess for ziter.
 * @return A pair indicating whether to continue iterating and the updated value of ziter.
 */
KOKKOS_FUNCTION
Kokkos::pair<bool, double> ImplicitIteration::iterate_rootfinding_algorithm(const double rprev,
                                                                            double ziter) const {
  // increment ziter
  const auto numerator = ode_gfunc(rprev, ziter);
  const auto denominator = ode_gfuncderivative(ziter);
  ziter = ziter * (1 - numerator / denominator);

  // test for next iteration
  const auto newnumerator = ode_gfunc(rprev, ziter);
  const auto do_iter = isnot_converged(newnumerator, numerator);

  return {do_iter, ziter};
}
