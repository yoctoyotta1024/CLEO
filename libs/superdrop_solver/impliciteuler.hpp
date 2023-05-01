// Author: Clara Bayley
// File: impliciteuler.hpp
/* Header for class implementing the implicit euler
  method for radial growth/shrink of each superdroplet
  due to condensation and diffusion of water
  vapour according to equations from "An Introduction
  To Clouds From The Microscale to Climate" by
  Lohmann, Luond and Mahrt, 1st edition." and
  Shima et al. 2009 */

#ifndef IMPLICITEULER_HPP
#define IMPLICITEULER_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "../claras_SDconstants.hpp"

namespace dlc = dimless_constants;

class ImplicitEuler
/* class for the Implicit Euler (IR) integration of
  to superdroplet condensational growth ODE. Uses a
  Newton Raphson (NR) root finding method to solve
  the implicit timestepping equation of stiff ODE */
{
private:
  const int maxiters; // maximum number of NR iterations before error raised
  const double delt;  // timestep of ODE solver (at each step NR iterations occur)
  const double rtol;  // relative tolerance for convergence of NR method
  const double atol;  // abolute tolerance for convergence of NR method

  struct IterReturn
  {
    bool do_iter;
    double ziter;
  };

  ImplicitEuler::IterReturn
  iterate_rootfinding_algorithm(double ziter, const double rprev,
                                const double s_ratio, const double akoh,
                                const double bkoh, const double ffactor) const;
  /* function performs one iteration of Newton Raphson rootfinding
  method and returns updated value of radius^2 alongside a boolean that
  is false if algorithm has converged */

  double initial_guess(const double s_ratio, const double akoh,
                      const double bkoh, const double r_k) const;
  /* returns appropriate initial value for ziter based on 
  uniqueness criteria of solution (root) of condensation ODE */

  double ode_gfunc(const double rsqrd, const double radius,
                   const double rprev, const double s_ratio,
                   const double akoh, const double bkoh,
                   const double ffactor) const;
  /* returns g(z) / z * delt for g(z) function used in root finding
  Newton Raphson Method for dr/dt condensation / evaporation ODE.
  ODE is for radial growth/shrink of each superdroplet due to
  condensation and diffusion of water vapour according to
  equations from "An Introduction To Clouds...."
  (see note at top of file). Note: z = ziter = radius^2 */

  double ode_gfuncderivative(const double rsqrd, const double radius,
                             const double s_ratio, const double akoh,
                             const double bkoh, const double ffactor) const;
  /* dg(z)/dz * delt, where dg(z)/dz is derivative of g(z) with
  respect to z=rsqrd. g(z) is polynomial to find root of using
  Newton Raphson Method. */

  inline bool isnotconverged(const double gfunciter,
                             const double gfuncprev) const
  /* boolean where True means
  criteria for ending newton raphson iteratiions
  has not yet been met. Criteria is standard local error
  test: |iteration - previous iteration|
  < RTOL * |iteration| + ATOL */
  {
    const double convergence_threshold = rtol * std::abs(gfunciter) + atol;
    const double currentvalue = std::abs(gfunciter - gfuncprev);

    return (currentvalue >= convergence_threshold); // true means it's not yet converged
  }

public:
  ImplicitEuler(const int maxiters, const double delt,
                const double rtol, const double atol)
      : maxiters(maxiters), delt(delt), rtol(rtol), atol(atol) {}

  double implicitmethod_forcondensation(const double s_ratio,
                                        const double akoh,
                                        const double bkoh,
                                        const double fkl,
                                        const double fdl,
                                        const double rprev) const;
  /* given initial guess for radius, (which is usually squared r from previous
  timestep 'rprev'^2), uses newton raphson iterative method to find value of r
  that converges on the root of function g(z), "gfunc", within the
  tolerances of the ImplicitEuler instance. Returns this value
  of r (usually used for new value of radius at current timestep)
  Refer to sect 5.1.2 Shima et al. 2009 for more details */
};

#endif // IMPLICITEULER_HPP
