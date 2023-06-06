// Author: Clara Bayley
// File: superdrop.cpp
/* Functionality file for superdrops class. Equations
	 referenced as (eqn [X.YY]) are from "An Introduction
	 To Clouds From The Microscale to Climate" by
	 Lohmann, Luond and Mahrt, 1st edition. */

#include "superdrop.hpp"

double Superdrop::mass() const
/* calculate total mass of droplet
	mass = (water + dry areosol)  */
{
	double mass(m_sol * (1.0 - solute->rho_l / solute->rho_sol)); // mass contribution of solute
	mass = 4.0 / 3.0 * M_PI * solute->rho_l * pow(radius, 3.0) + mass;

	return mass;
}

double Superdrop::equilibrium_wetradius(const double s_ratio,
                                        const double temp) const
/* Performs Newton Raphson root finding algorithm using functions in 
WetRadius root finder struct to solve equation for equilibrium (wet) 
radius of superdroplet at given relative humidity. Equilibrium radius 
defined by radius when ODE from "An Introduction To Clouds...."
(see note at top of file) eqn [7.28] = 0. */
{
  const unsigned int maxiters = 100;
  const double akoh = akohler_factor(temp);
  const double bkoh = bkohler_factor();

	WetRadius wrrf{maxiters};

  return wrrf.get_wetradius(radius, s_ratio, akoh, bkoh);
}

double Superdrop::rhoeff() const
/* calculates effective density of droplet
so mass_droplet = m = 4/3*pi*r^3 * rhoeff */
{
	double effsol(1.0 - solute->rho_l / solute->rho_sol);
	effsol = 3.0 * m_sol / (4.0 * M_PI * pow(radius, 3.0)) * effsol; // effect of solute on density

	return solute->rho_l + effsol;
}

// double Superdrop::mass_liq() const
// /* mass of only water in droplet */
// {
// 	return solute->rho_l * vol_liq();
// }

double Superdrop::akohler_factor(const double temp) const
/* calculate value of a in raoult factor (exp^(a/r))
	to account for effect of dissolved solute
	on radial growth of droplet. Using equations from
	"An Introduction To Clouds...." (see note at top of file) */
{
	constexpr double akoh = 3.3e-7 / (dlc::TEMP0 * dlc::R0);

	return akoh / temp; // dimensionless version of eqn [6.24]
}

double Superdrop::bkohler_factor() const
/* calculate value of b in kelvin factor (1-b/r^3)
	to account for curvature on radial growth
	of droplet. Using equations from "An Introduction
	To Clouds...." (see note at top of file) */
{
	constexpr double bkoh = 4.3e-6 * dlc::RHO0 / dlc::MR0;

	return bkoh * m_sol * solute->ionic / solute->mrsol; // dimensionless version of eqn [6.22]
}

double Superdrop::change_radius(const double newradius)
/* Update droplet radius to newradius or dry_radius() and
return resultant change in radius (delta_radius = newradius-radius). 
Prevents drops shrinking further once they are size of dry_radius(). */
{
	/*  if droplets are dry, do not shrink further */
	const double oldradius = radius;
	radius = std::max(dry_radius(), newradius);

	/* return change in radius due to growth/shrinking of droplet */
	return radius - oldradius;
}

double WetRadius::get_wetradius(const double radius0, const double s_ratio,
                                const double akoh, const double bkoh) const
/* Iterate Newton Raphson root finding algorithm to
return wet radius of a superdroplet in equilibrium
with supersaturation s_ratio */
{
  unsigned int iter = 0;
  bool do_iter = true;
  double ziter = radius0; // value of ziter at iter=0 (no iterations yet)
  
  while (do_iter)
  {
    if (iter > maxiters)
    {
      const std::string err("Newton Raphson Method did not converge"
                            " within " +
                            std::to_string(maxiters) +
                            " iterations to find wet radius\n");
      throw std::invalid_argument(err);             
    }
    else
    {
      /* add one to the number of attempted iterations
        z^(m+1) for iteration m+1 starting at m=0 */
      iter += 1;
      const auto a(iterate_rootfinding(ziter, s_ratio, akoh, bkoh));
      do_iter = a.do_iter;
      ziter = a.ziter;
    }
  }

  return ziter;
}

WetRadius::IterReturn
WetRadius::iterate_rootfinding(double ziter, const double s_ratio,
                               const double akoh, const double bkoh) const
/* performs 1 iteration of newton raphson root finding algorithm for 
obtaining the equilibrium wet radius of the condensation ODE at a given 
relative humidity (s_ratio). ODE from "An Introduction To Clouds...." 
(see note at top of file) eqn [7.28] */
{
  const double ode = wetradius_polynomial(ziter, s_ratio, akoh, bkoh);
  const double odederiv = 3 * (s_ratio - 1.0) * pow(ziter, 2.0) - 2 * akoh * ziter;
  
  ziter -= ode / odederiv; // increment ziter

  // prepare for next iteration or end while loop
  const double new_ode = wetradius_polynomial(ziter, s_ratio, akoh, bkoh); 
  
  return IterReturn{isnotconverged(new_ode, ode), ziter};
}

double WetRadius::wetradius_polynomial(const double ziter,
                                       const double s_ratio,
                                       const double akoh,
                                       const double bkoh) const
/* returns value of (cubic) polynomial evaluted at ziter. Root of this
polynomial is the value of the equilibrium (wet) radius of a superdroplet
at a given relative humidity (ie. s_ratio) derived from ODE in
"An Introduction To Clouds...." (see note at top of file) eqn [7.28] */
{
  return (s_ratio - 1.0) * pow(ziter, 3.0) - (akoh * pow(ziter, 2.0)) + bkoh;
}

bool WetRadius::isnotconverged(const double new_ode,
                               const double ode) const
/* boolean where True means criteria for ending newton raphson iterations
has not yet been met. Criteria is standard local error test:
|iteration - previous iteration| < RTOL * |iteration| + ATOL */
{
  const double rtol = 1e-8;
  const double atol = 1e-8;

  const double convergence_threshold = rtol * std::abs(new_ode) + atol;
  const double currentvalue = std::abs(new_ode - ode);

  return (currentvalue >= convergence_threshold); // true means it's not yet converged
}