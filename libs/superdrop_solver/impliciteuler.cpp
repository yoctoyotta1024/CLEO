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

double ImplicitEuler::solve_condensation(const double s_ratio,
                                         const double akoh,
                                         const double bkoh,
                                         const double fkl,
                                         const double fdl,
                                         const double rprev) const
/* forward timestep previous radius 'rprev' by delt using an implicit
euler method to integrate the condensation/evaporation ODg. Implict
timestepping equation defined in section 5.1.2 of Shima et al. 2009
and is root of polynomial g(z) = 0, where z = [R_i(t+delt)]^squared.
Newton Raphson iterations are used to converge towards the root of 
g(z) within the tolerances of an ImpIter instance. Tolerances,
maxium number of iterations and sub-timestepping are adjusted based on
on the uniqueness criteria of the polynomial g(z). Refer to section
5.1.2 Shima et al. 2009 and section 3.3.3 of Matsushima et al. 2023
for more details. */
{
  const double ffactor(dlc::Rho_l * (fkl + fdl));
  const double ract_ratio(rprev * rprev * akoh / 3.0 / bkoh);

  const double max_uniquedelt(
      12.5 * ffactor / akoh * std::pow(5.0 * bkoh / akoh, 1.5));

// std::cout << "-----\n";
  if ((s_ratio <= 1.0 && ract_ratio < 1.0) || (delt <= max_uniquedelt))
  {
    return unique_solution_for_implicitmethod(s_ratio, akoh, bkoh,
                                              ffactor, rprev);
  }
  else
  {
    return subtimestep_solution_for_implicitmethod(s_ratio, akoh, bkoh,
                                                   ffactor, rprev);
  }
}

double ImplicitEuler::initial_ziter(const double rprev,
                                    const double s_ratio,
                                    const double akoh,
                                    const double bkoh) const
/* returns appropriate initial value (ie. a reasonable guess) for
'ziter' to use as first iteration of newton raphson method in
rootfinding algorithm for timestepping condensation/evaporation ODE */
{
  // const double s_activ(1 + std::sqrt(4.0*std::pow(akoh, 3.0) / (27*bkoh))); // activation supersaturation of droplet
  // if (s_ratio > s_activ)
  // {
  //   return std::pow(1e-3 / dlc::R0, 2.0);
  // }
  
  // const double r1sqrd(bkoh/akoh); // (equilibrium radius for drolet at s_ratio=1)^2
  
  // return std::max(std::pow(rprev, 2.0), r1sqrd);
  return std::pow(rprev, 2.0);
}

double ImplicitEuler::
    unique_solution_for_implicitmethod(const double s_ratio,
                                       const double akoh,
                                       const double bkoh,
                                       const double ffactor,
                                       const double rprev) const
/* construct ImpIter instance to timestep condensation ODE
by delt assuming that solution to g(ziter)=0 is unique and therefore
Newton Raphson root finding algorithm converges quickly. This
means sufficiently small tolerances and timestep are comparitively
large, and maximum number of iterations is small:
Relative tolerance 'rtol' >= 0.1, absolute tolerance 'atol' >= 0.1.
Maximum number of Newton Raphson Iterations for timestep delt <= 10 */
{
  const unsigned int iterlimit(std::min(10, maxiters)); // at most 10 iterations
  const double subdelt(delt);                  // no subtimestepping
  const double rtol(std::max(1.0, maxrtol));   // at least 0.1 rtol
  const double atol(std::max(1.0, maxatol));   // at least 0.1 atol

  const ImpIter impit{iterlimit, subdelt, rtol, atol,
                      s_ratio, akoh, bkoh, ffactor, rprev};
  
  const double init_ziter(std::pow(rprev, 2.0)); // guess for ziter (before any iterations)
  return impit.implicitmethod_forcondensation(init_ziter);
}

double ImplicitEuler::
    subtimestep_solution_for_implicitmethod(const double s_ratio,
                                            const double akoh,
                                            const double bkoh,
                                            const double ffactor,
                                            const double rprev) const
/* construct ImpIter instance to timestep condensation ODE
by delt assuming that solution to g(ziter)=0 is not unique and
therefore sub-timestepping is required with sufficiently
small tolerances. For each subtimestep, perform Newton Raphson
root finding algorithm with comparitevly small tolerances to
obtain solution to g(ziter) polynomial. Tolerances are:
relative tolerance 'rtol' <= 0.01, absolute tolerance 'atol' <= 0.1.
Subtimestep = delt/10 and maximum number of Newton Raphson
Iterations >= 25 for each subtimetep. */
{
  const int iterlimit(std::max(25, maxiters));
  const double subdelt(delt);
  const double rtol(std::min(0.01, rtol));
  const double atol(std::min(0.1, atol));

  const ImpIter impit{iterlimit, subdelt, rtol, atol,
                      s_ratio, akoh, bkoh, ffactor, rprev};

  const double init_ziter(initial_ziter(rprev, s_ratio,
                                        akoh, bkoh)); // guess for ziter (before any iterations)
  return impit.implicitmethod_forcondensation(init_ziter);
}

double ImplicitEuler::ImpIter::
    implicitmethod_forcondensation(double ziter) const
/* given initial guess for ziter, (which is usually radius^squared
from previous timestep), uses newton raphson iterative method to
find new value of radius that converges on the root of the
polynomial g(ziter) within the tolerances of the ImpIter instance.
Raises error if method does not converge within 'iterlimit' iterations.
Otherwise returns new value for the radius (which is the radius at
timestep 't+subdelt'. Refer to section 5.1.2 Shima et al. 2009
and section 3.3.3 of Matsushima et al. 2023 for more details. */
{
  bool do_iter(true);
  unsigned int iter(1);

  while (do_iter)
  {
    if (iter > iterlimit)
    {
      const std::string err = "Newton Raphson Method did not converge "
                          "within " + std::to_string(iterlimit) + " iterations\n";
      throw std::invalid_argument(err);
      break;
    }
    else
    {
      /* perform one attempted iteration  ziter^(m) -> ziter^(m+1)
      for iteration m+1 starting at m=1 */
      const auto iterret = iterate_rootfinding_algorithm(ziter);
      do_iter = iterret.first;
      ziter = iterret.second;
      iter += 1;
    }
  }

  // once iterations have converged, return ziter
  return std::sqrt(ziter);
}

std::pair<bool, double> ImplicitEuler::ImpIter::
    iterate_rootfinding_algorithm(double ziter) const
/* function performs one iteration of Newton Raphson rootfinding
  method and returns updated value of radius^2 alongside a boolean that
  is false if algorithm has converged */
{ 
  // increment ziter
  const double numerator = ode_gfunc(ziter);
  const double denominator = ode_gfuncderivative(ziter);
  ziter = ziter * (1 - numerator / denominator); 

  // test for next iteration
  const double newnumerator = ode_gfunc(ziter);
  const bool do_iter = isnotconverged(newnumerator, numerator);

  return std::pair<bool, double>{do_iter, ziter};
}

double ImplicitEuler::ImpIter::
    ode_gfunc(const double rsqrd) const
/* returns g(z) / z * delt for g(z) function used in root finding
Newton Raphson Method for dr/dt condensation / evaporation ODE. 
ODE is for radial growth/shrink of each superdroplet due to
condensation and diffusion of water vapour according to
equations from "An Introduction To Clouds...."
(see note at top of file). Note: z = ziter = radius^2 */
{
  const double radius = std::sqrt(rsqrd);

  const double alpha(s_ratio - 1 - akoh / radius + bkoh / std::pow(radius, 3.0));
  const double beta(2.0 * subdelt / (rsqrd * ffactor));
  const double gamma(std::pow(rprev/radius, 2.0));

  return 1 - gamma - alpha * beta;
}

double ImplicitEuler::ImpIter::
    ode_gfuncderivative(const double rsqrd) const
/* dg(z)/dz * delt, where dg(z)/dz is derivative of g(z) with
respect to z=rsqrd. g(z) is polynomial to find root of using
Newton Raphson Method. */
{
  const double radius = std::sqrt(rsqrd);
  
  const double alpha(akoh/radius - 3.0 * bkoh/ std::pow(radius, 3.0));
  const double beta(subdelt / (rsqrd * ffactor));

  return 1 - alpha * beta;
}