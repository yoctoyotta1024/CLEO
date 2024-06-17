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
 * Last Modified: Monday 17th June 2024
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
 * @brief Struct for performing iterations of the Implicit Euler Method.
 *
 * This struct defines parameters and methods for performing an iterations of the implicit method.
 *
 * _Note: abbreviation NR = Newton Raphson (Method)
 *
 */
struct ImplicitIterations {
 private:
  size_t maxniters; /**< Maximum no. iterations of Newton Raphson Method */
  double rtol;      /**< Relative tolerance for convergence of NR method. */
  double atol;      /**< Absolute tolerance for convergence of NR method. */

 public:
  /**
   * @brief Struct for constants of ODE during integration
   */
  struct ODEConstants {
    double s_ratio;  ///< Supersaturation ratio.
    double akoh;     ///< Kelvin factor in Kohler theory "a".
    double bkoh;     ///< Raoult factor in Kohler theory "b".
    double ffactor;  ///< Sum of heat and vapor diffusion factors.
  };

  ImplicitIterations(const size_t maxniters, const double rtol, const double atol)
      : maxniters(maxniters), rtol(rtol), atol(atol) {}

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
   * @param subdelt Time over which to integrate ODE
   * @param rprev Radius of droplet at previous timestep.
   */
  KOKKOS_FUNCTION double integrate_condensation_ode(const ODEConstants &odeconsts,
                                                    const double subdelt, const double rprev) const;

 private:
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
  KOKKOS_FUNCTION double initialguess(const ODEConstants &odeconsts, const double rprev) const;
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
  KOKKOS_FUNCTION Kokkos::pair<double, bool> newtonraphson_niterations(
      const ODEConstants &odeconsts, const double subdelt, const double rprev, double ziter,
      const size_t niters) const;
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
   * @param maxniters The maxiumum number of iterations to attempt.
   * @param rprev Radius at the previous timestep.
   * @param ziter The current guess for ziter.
   * @return The updated value of ziter.
   */
  KOKKOS_FUNCTION double newtonraphson_untilconverged(const ODEConstants &odeconsts,
                                                      const size_t maxniters, const double subdelt,
                                                      const double rprev, double ziter) const;

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
  Kokkos::pair<double, bool> iterate_rootfinding_algorithm(const ODEConstants &odeconsts,
                                                           const double subdelt, const double rprev,
                                                           double ziter) const;

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
  KOKKOS_FUNCTION double ode_gfunc(const ODEConstants &odeconsts, const double subdelt,
                                   const double rprev, const double rsqrd) const;

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
  KOKKOS_FUNCTION double ode_gfuncderivative(const ODEConstants &odeconsts, const double subdelt,
                                             const double rsqrd) const;

  /**
   * @brief Returns true if the Newton-Raphson iterations have converged.
   *
   * This method checks if the Newton-Raphson iterations have converged based on a standard
   * local error test: |iteration - previous iteration| < RTOL * |iteration| + ATOL.
   *
   * @param gfunciter value proportional to g(z) for current iteration
   * @param gfuncprev value proportional to g(z) for previous iteration
   * @return boolean=true if converged, false otherwise
   */
  KOKKOS_FUNCTION
  bool check_for_convergence(const double gfunciter, const double gfuncprev) const {
    const auto converged = double{rtol * Kokkos::abs(gfunciter) + atol};
    const auto currentvalue = double{Kokkos::abs(gfunciter - gfuncprev)};

    return (currentvalue < converged);  // true means converged
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
  double delt;       /**< Timestep of ODE solver (at each step implicit method is called). */
  double minsubdelt; /**< Minimum subtimestep in cases of substepping */
  ImplicitIterations implit; /**< Performs Newton Raphson Iterations of Implicit Method */

 public:
  /**
   * @brief Constructor for ImplicitEuler class.
   */
  ImplicitEuler(const size_t maxniters, const double delt, const double rtol, const double atol,
                const double minsubdelt)
      : delt(delt), minsubdelt(minsubdelt), implit(maxniters, rtol, atol) {
    // TODO(CB): WIP new config params to match this constructor and re-order constructor args
  }

  /**
   * @brief Integrates the condensation / evaporation ODE employing the Implicit Euler method
   * similarly to Matsushima et. al, 2023.
   *
   * Forward timestep previous radius 'rprev' by delt using an Implicit Euler method (possibly with
   * sub-timestepping) to integrate the condensation/evaporation ODE using fixed thermodynamics from
   * the start of the timestep.
   *
   * @param s_ratio The saturation ratio.
   * @param kohler_ab A pair containing 'a' and 'b' factors for Kohler curve in that order.
   * @param ffactor The sum of the diffusion factors.
   * @param rprev Previous radius at time = t
   * @return Updated radius for time = t + delt
   */
  KOKKOS_FUNCTION double solve_condensation(const double s_ratio,
                                            const Kokkos::pair<double, double> kohler_ab,
                                            const double ffactor, const double rprev) const;
};

#endif  // LIBS_SUPERDROPS_IMPLICITEULER_HPP_
