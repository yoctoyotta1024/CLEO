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
  iterate_rootfinding_algorithm(double ziter, const double s_ratio,
                                const double akoh, const double bkoh,
                                const double fkl, const double fdl,
                                const double r_k) const;
  /* function performs one iteration of Newton Raphson rootfinding
  method and returns updated value of radius alongside a boolean that
  is false if algorithm has converged */

  double ode_gfunc(const double ziter, const double r_k,
                   const double s_ratio, const double akoh,
                   const double bkoh, const double fkl,
                   const double fdl) const;
  /* returns value of g(z) function
  used in root finding Newton Raphson
  Method for condensation ODE */

  double condensation_ode(const double radius, const double s_ratio,
                          const double akoh, const double bkoh,
                          const double fkl, const double fdl) const;
  /* dr/dt ODE for radial growth/shrink
   of each superdroplet due to	condensation and
   diffusion of water vapour according to
   equations from "An Introduction To
   Clouds...." (see note at top of file)
   used in calculation of ode_gfunc(...) for
   Newton Raphson method increment */

  double ode_gfuncderivative(const double ziter, const double s_ratio,
                             const double akoh, const double bkoh,
                             const double fkl, const double fdl) const;
  /* derivative of g(z) function w.r.t used
  in root finding Newton Raphson Method.
  Involves calculation of r_rdot_deriv,
  which is deriative of r * dr/dt ODE */

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

  double newton_raphson_iterator_forcondensation(const double s_ratio,
                                                 const double akoh, const double bkoh,
                                                 const double fkl, const double fdl,
                                                 const double r_k) const;
  /* given initial guess for radius, "r_k", (which is usually r from
  previous timestep), newton_raphson_iterator finds value of r
  that converges on the root of function g(z), "gfunc", within the
  tolerances of the ImplicitEuler instance. Returns this value
  of r (usually used for new value of radius at current timestep).
  Refer to sect 5.1.2 Shima et al. 2009 for more details */
};

#endif // IMPLICITEULER_HPP
