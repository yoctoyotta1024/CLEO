/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: impliciteuler.hpp
 * Project: superdrops
 * Created Date: Thursday 26th October 2023
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
 * Header for class implementing the Implicit Euler method for radial growth/shrink of each
 * droplet due to condensation / evaporation and diffusion of water vapour according
 * to equations from "An Introduction To Clouds From The Microscale to Climate" by Lohmann,
 * Luond and Mahrt, 1st edition." and Shima et al. 2009
 */

#ifndef LIBS_SUPERDROPS_IMPLICITEULER_HPP_
#define LIBS_SUPERDROPS_IMPLICITEULER_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <cassert>

#include "../cleoconstants.hpp"

namespace dlc = dimless_constants;

/**
 * @brief Struct defining parameters and methods for iterations of the Implicit Euler Method.
 *
 * This struct defines parameters and methods for performing an iterations of the implicit method.
 * E.g. struct defines the number of iterations to try before testing convergence, the timestep of
 * the method, it's relative and absolute tolerances, and methods for computing the initial guess
 * for and performing the Newton-Raphson root finding.
 *
 * _Note: abbreviation NR = Newton Raphson (Method)
 *
 */
struct ImplicitIteration {
  unsigned int niters;  ///< Number of NR iterations to try before testing convergence.
  double subdelt;  ///< Timestep of implicit method (at each substep >= niter NR iterations occur).
  double rtol;     ///< Relative tolerance for convergence of NR method.
  double atol;     ///< Absolute tolerance for convergence of NR method.
  double s_ratio;  ///< Supersaturation ratio.
  double akoh;     ///< Kelvin factor in Kohler theory "a".
  double bkoh;     ///< Raoult factor in Kohler theory "b".
  double ffactor;  ///< Sum of heat and vapor diffusion factors.

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
  KOKKOS_FUNCTION double initialguess(const double rprev) const;

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
  KOKKOS_FUNCTION
  double initialguess_shima(const double rprev) const;

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
  KOKKOS_FUNCTION
  double newtonraphson_niterations(const double rprev, double ziter) const;

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
  KOKKOS_FUNCTION
  double newtonraphson_untilconverged(const unsigned int iterlimit, const double rprev,
                                      double ziter) const;

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
  KOKKOS_FUNCTION Kokkos::pair<bool, double> iterate_rootfinding_algorithm(const double rprev,
                                                                           double ziter) const;

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
  KOKKOS_FUNCTION
  double ode_gfunc(const double rprev, const double rsqrd) const;

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
  KOKKOS_FUNCTION
  double ode_gfuncderivative(const double rsqrd) const;

  /**
   * @brief Returns true if the Newton-Raphson iterations have *not* yet converged.
   *
   * This method checks if the Newton-Raphson iterations have converged based on a standard
   * local error test: |iteration - previous iteration| < RTOL * |iteration| + ATOL.
   *
   * @param gfunciter g(z) of current iteration
   * @param gfuncprev  g(z) of previous iteration
   * @return boolean=true if not yet converged, false otherwise
   */
  KOKKOS_FUNCTION
  bool isnot_converged(const double gfunciter, const double gfuncprev) const {
    const auto converged = double{rtol * Kokkos::abs(gfunciter) + atol};
    const auto currentvalue = double{Kokkos::abs(gfunciter - gfuncprev)};

    return (currentvalue >= converged);  // true means it's not yet converged
  }
};

/**
 * @brief Class for Implicit Euler (IE) integration of superdroplet condensational growth /
 * evaporational shrinking ODE.
 *
 * This class performs Implicit Euler integration of the superdroplet condensation / evaporation ODE
 * using a Newton-Raphson root finding method to solve the implicit timestep equation of a stiff
 * ODE.
 */
class ImplicitEuler {
 private:
  unsigned int niters;
  /**< Suggested number of iterations for Newton Raphson method before testing for convergence. */
  double delt;    /**< Timestep of ODE solver (at each step implicit method is called). */
  double maxrtol; /**< Adjustable relative tolerance for convergence of NR method. */
  double maxatol; /**< Adjustable absolute tolerance for convergence of NR method. */
  double subdelt; /**< Number of substeps to intergation when supersat close to 1 (0 = false). */

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
  KOKKOS_FUNCTION
  double substepped_implicitmethod(const ImplicitIteration &implit, const unsigned int nsubsteps,
                                   const double rprev) const;

 public:
  /**
   * @brief Constructor for ImplicitEuler class.
   */
  ImplicitEuler(const unsigned int niters, const double delt, const double maxrtol,
                const double maxatol, const double subdelt)
      : niters(niters), delt(delt), maxrtol(maxrtol), maxatol(maxatol), subdelt(subdelt) {}

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
  KOKKOS_FUNCTION
  double solve_condensation_matsushima(const double s_ratio,
                                       const Kokkos::pair<double, double> akoh_bkoh,
                                       const double ffactor, const double rprev) const;

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
  KOKKOS_FUNCTION
  double solve_condensation(const double s_ratio, const Kokkos::pair<double, double> akoh_bkoh,
                            const double ffactor, const double rprev) const;
};

#endif  // LIBS_SUPERDROPS_IMPLICITEULER_HPP_
