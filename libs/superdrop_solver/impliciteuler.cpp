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
  const double ffactor(dlc::Rho_l * (fkl + fdl));

  double ziter(initial_guess(rprev, s_ratio, akoh, bkoh)); // ziter at iter=0 (before any iterations)
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
      const IterReturn a = iterate_rootfinding_algorithm(ziter, rprev, s_ratio,
                                                         akoh, bkoh, ffactor);
      do_iter = a.do_iter;
      ziter = a.ziter;
      iter += 1;
    }
  }

  // once iterations have converged, return ziter
  return std::sqrt(ziter);
}

double ImplicitEuler::initial_guess(const double rprev,
                                    const double s_ratio,
                                    const double akoh,
                                    const double bkoh) const
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
                                             const double rprev,
                                             const double s_ratio,
                                             const double akoh,
                                             const double bkoh,
                                             const double ffactor) const
/* function performs one iteration of Newton Raphson rootfinding
  method and returns updated value of radius^2 alongside a boolean that
  is false if algorithm has converged */
{ 
  const double radius = std::sqrt(ziter);

  // increment ziter
  const double numerator = ode_gfunc(ziter, radius, rprev, s_ratio,
                                     akoh, bkoh, ffactor);
  const double denominator = ode_gfuncderivative(ziter, radius, s_ratio,
                                                 akoh, bkoh, ffactor);
  ziter = ziter * (1 - numerator / denominator); 

  // test for next iteration
  const double newnumerator = ode_gfunc(ziter, radius, rprev, s_ratio,
                                        akoh, bkoh, ffactor);
  const bool do_iter = isnotconverged(newnumerator, numerator);

  return IterReturn{do_iter, ziter};
}

double ImplicitEuler::ode_gfunc(const double rsqrd, const double radius,
                                const double rprev, const double s_ratio,
                                const double akoh, const double bkoh,
                                const double ffactor) const
/* returns g(z) / z * delt for g(z) function used in root finding
Newton Raphson Method for dr/dt condensation / evaporation ODE. 
ODE is for radial growth/shrink of each superdroplet due to
condensation and diffusion of water vapour according to
equations from "An Introduction To Clouds...."
(see note at top of file). Note: z = ziter = radius^2 */
{
  const double alpha(s_ratio - 1 - akoh / radius + bkoh / std::pow(radius, 3.0));

  const double beta(2.0 * delt / (rsqrd * ffactor));

  const double gamma(std::pow(rprev/radius, 2.0));

  return 1 - gamma - alpha * beta;
}

double ImplicitEuler::ode_gfuncderivative(const double rsqrd,
                                          const double radius,
                                          const double s_ratio,
                                          const double akoh,
                                          const double bkoh,
                                          const double ffactor) const
/* dg(z)/dz * delt, where dg(z)/dz is derivative of g(z) with
respect to z=rsqrd. g(z) is polynomial to find root of using
Newton Raphson Method. */
{
  const double alpha(akoh/radius - 3.0 * bkoh/ std::pow(radius, 3.0));

  const double beta(delt / (rsqrd * ffactor));

  return 1 - alpha * beta;
}