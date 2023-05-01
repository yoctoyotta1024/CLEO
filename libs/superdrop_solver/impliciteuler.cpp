// Author: Clara Bayley
// File: impliciteuler.cpp
/* Class implementing the implicit euler
  method for radial growth/shrink of each superdroplet
  due to condensation and diffusion of water
  vapour according to equations from "An Introduction
  To Clouds From The Microscale to Climate" by
  Lohmann, Luond and Mahrt, 1st edition." and
  Shima et al. 2009 */

#include "impliciteuler.hpp"

double ImplicitEuler::implicitmethod_forcondensation(const double s_ratio,
                                                     const double akoh,
                                                     const double bkoh,
                                                     const double fkl,
                                                     const double fdl,
                                                     const double rprev) const
/* given initial guess for radius, (which is usually squared r from previous
timestep 'rprev'^2), uses newton raphson iterative method to find value of r
that converges on the root of function g(z), "gfunc", within the
tolerances of the ImplicitEuler instance. Returns this value
of r (usually used for new value of radius at current timestep)
Refer to sect 5.1.2 Shima et al. 2009 for more details */
{
  double ziter(initial_guess(s_ratio, akoh, bkoh, rprev)); // ziter at iter=0 (before any iterations)
  bool do_iter(true);
  int iter(1);

  while (do_iter)
  {
    if (iter > maxiters)
    {
      const std::string err = "Newton Raphson Method did not converge "
                          "within " + std::to_string(maxiters) + " iterations\n";
      throw std::invalid_argument(err);
      break;
    }
    else
    {
      /* add one to the number of attempted iterations
        z^(m+1) for iteration m+1 starting at m=0 */
      const IterReturn a = iterate_rootfinding_algorithm(ziter, s_ratio, akoh,
                                                         bkoh, fkl, fdl, rprev);
      do_iter = a.do_iter;
      ziter = a.ziter;
      iter += 1;
    }
  }

  // once iterations have converged, return ziter
  return std::sqrt(ziter);
}

double ImplicitEuler::initial_guess(const double s_ratio, const double akoh,
                                    const double bkoh, const double rprev) const
/* returns appropriate initial value for ziter based on 
uniqueness criteria of solution (root) of condensation ODE */
{
  // const double sact_factor = 4.0*std::pow(akoh, 3.0) / (27*bkoh);
  // if (s_ratio > 1.0 && sact_factor < 1.0)
  // {
  //   return std::pow(1e-6 / dlc::R0, 2.0);
  // }

  const double r1sqrd(bkoh/akoh); // (equilibrium radius for drolet at s_ratio=1)^2
  
  return std::max(std::pow(rprev, 2.0), r1sqrd);
}

ImplicitEuler::IterReturn
ImplicitEuler::iterate_rootfinding_algorithm(double ziter,
                                             const double s_ratio,
                                             const double akoh,
                                             const double bkoh,
                                             const double fkl,
                                             const double fdl,
                                             const double rprev) const
/* function performs one iteration of Newton Raphson rootfinding
  method and returns updated value of radius^2 alongside a boolean that
  is false if algorithm has converged */
{ 
  const double radius = std::sqrt(ziter);

  const double numerator = ode_gfunc(radius, rprev, s_ratio, akoh, bkoh, fkl, fdl);
  const double denominator = ode_gfuncderivative(radius, s_ratio, akoh, bkoh, fkl, fdl);
  ziter -= ziter * numerator / denominator; // increment ziter

  // prepare for next iteration or end while loop
  const double newnumerator = ode_gfunc(radius, rprev, s_ratio, akoh, bkoh, fkl, fdl);
  const bool do_iter = isnotconverged(newnumerator, numerator);

  return IterReturn{do_iter, ziter};
}

double ImplicitEuler::ode_gfunc(const double ziter, const double r_k,
                                const double s_ratio, const double akoh,
                                const double bkoh, const double fkl,
                                const double fdl) const
/* returns value of g(z) function
used in root finding Newton Raphson
Method for condensation ODE */
{
  const double rdot = condensation_ode(ziter, s_ratio, akoh, bkoh, fkl, fdl);

  const double gfunc = (pow(ziter, 2.0) - pow(r_k, 2.0)) / (2 * delt) - ziter * rdot;

  return gfunc;
}

double ImplicitEuler::condensation_ode(const double radius,
                                       const double s_ratio,
                                       const double akoh,
                                       const double bkoh,
                                       const double fkl,
                                       const double fdl) const
/* dr/dt ODE for radial growth/shrink
  of each superdroplet due to	condensation and
  diffusion of water vapour according to
  equations from "An Introduction To
  Clouds...." (see note at top of file)
  used in calculation of ode_gfunc(...) for
  Newton Raphson method increment */
{
  const double nominator = (s_ratio - 1.0) - (akoh / radius) + (bkoh / pow(radius, 3.0)); // eqn [7.28]

  const double denominator = dlc::Rho_l * (fkl + fdl) * radius;

  const double rdot = nominator / denominator;

  return rdot;
}

double ImplicitEuler::ode_gfuncderivative(const double ziter,
                                          const double s_ratio,
                                          const double akoh,
                                          const double bkoh,
                                          const double fkl,
                                          const double fdl) const
/* derivative of g(z) function w.r.t z.
used in root finding Newton Raphson Method.
Involves calculation of r_rdot_deriv,
which is deriative of r * dr/dt ODE */
{
  const double numerator = (akoh / pow(ziter, 2.0)) - (3 * bkoh / pow(ziter, 4.0));
  const double denominator = dlc::Rho_l * (fkl + fdl);
  const double r_rdot_deriv = numerator / denominator; // derivative of r * eqn [7.28] (make sure this matches with condensation_ode)

  const double gfuncderiv = ziter / delt - r_rdot_deriv;

  return gfuncderiv;
}