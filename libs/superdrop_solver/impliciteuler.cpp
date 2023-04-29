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

double ImplicitEuler::newton_raphson_iterator_forcondensation(const double s_ratio,
                                                              const double akoh,
                                                              const double bkoh,
                                                              const double fkl,
                                                              const double fdl,
                                                              const double r_k) const
/* given initial guess for radius, "r_k", (which is usually r from
previous timestep), newton_raphson_iterator finds value of r
that converges on the root of function g(z), "gfunc", within the
tolerances of the ImplicitEuler instance. Returns this value
of r (usually used for new value of radius at current timestep)
Refer to sect 5.1.2 Shima et al. 2009 for more details */
{
  bool do_iter = true;
  int iter = 0;
  double ziter(initial_guess(s_ratio, akoh, bkoh, r_k)); // ziter at iter=0 (before any iterations)

  while (do_iter)
  {
    if (iter > maxiters)
    {
      const std::string errormsg = "WARNING! Newton Raphson Method did not "
                          " converge within " +
                          std::to_string(maxiters) +
                          " no. iterations\n";
      throw std::invalid_argument(errormsg);
      break;
    }
    else
    {
      /* add one to the number of attempted iterations
        z^(m+1) for iteration m+1 starting at m=0 */
      iter += 1;
      const IterReturn a = iterate_rootfinding_algorithm(ziter, s_ratio, akoh,
                                                         bkoh, fkl, fdl, r_k);
      do_iter = a.do_iter;
      ziter = a.ziter;
    }
  }

  // once iterations have converged, return ziter
  return ziter;
}

double ImplicitEuler::initial_guess(const double s_ratio, const double akoh,
                                    const double bkoh, const double r_k) const
/* returns appropriate initial value for ziter based on 
uniqueness criteria of solution (root) of condensation ODE */
{
  // const double sact_factor = 4.0*std::pow(akoh, 3.0) / (27*bkoh);
  // if (s_ratio > 1.0 && sact_factor < 1.0)
  // {
  //   return 1e-6 / dlc::R0;
  // }

  const double r1(std::sqrt(bkoh/akoh)); // equilibrium radius for drolet at s_ratio=1
  
  return std::max(r_k, r1);
}

ImplicitEuler::IterReturn
ImplicitEuler::iterate_rootfinding_algorithm(double ziter,
                                             const double s_ratio,
                                             const double akoh,
                                             const double bkoh,
                                             const double fkl,
                                             const double fdl,
                                             const double r_k) const
/* function performs one iteration of Newton Raphson rootfinding
  method and returns updated value of radius alongside a boolean that
  is false if algorithm has converged */
{ 
  const double gfunc = ode_gfunc(ziter, r_k, s_ratio, akoh, bkoh, fkl, fdl);
  const double gfuncderiv = ode_gfuncderivative(ziter, s_ratio, akoh, bkoh, fkl, fdl);
  ziter -= gfunc / gfuncderiv; // increment ziter

  // prepare for next iteration or end while loop
  const double newgfunc = ode_gfunc(ziter, r_k, s_ratio, akoh, bkoh, fkl, fdl);
  const bool do_iter = isnotconverged(newgfunc, gfunc);

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