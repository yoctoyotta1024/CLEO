// Author: Clara Bayley
// File: run_cvodesdm.cpp
/* Functions involved specifically
in running SDM coupled to Sundials
CVODE ODE solver for the thermodynamics.
Coupling can be one-way or both ways */

#include "./run_cvodesdm.hpp"

std::vector<double> initcvodethermo(const size_t num_gridboxes,
                                        const Config &config)
/* return vector of dimensionless initial conditions
for thermodyanmic variables (p, temp, qv, qc) to
initialise cvode thermodynamics solver */
{
  constexpr int NVARS = 4;                  // no. (distinct) variables per grid box
  const size_t neq = NVARS * num_gridboxes; // total no. variables
  std::vector<double> y_init(neq);

  const double p_init = config.P_INIT / dlc::P0;
  const double temp_init = config.TEMP_INIT / dlc::TEMP0;
  const double vapourp_init = saturation_pressure(temp_init) * config.relh_init / 100.0;
  const double qv_init = vapourpressure_2_massmixratio(vapourp_init, p_init);
  const double qc_init = config.qc_init;

  for (size_t k = 0; k < neq; k += NVARS)
  {
    y_init[k] = p_init;
    y_init[k + 1] = temp_init;
    y_init[k + 2] = qv_init;
    y_init[k + 3] = qc_init;
  }

  return y_init;
}

std::mt19937 preparetotimestep(CvodeThermoSolver &cvode,
                               std::vector<GridBox> &gridboxes,
                               const bool wetradiiinit,
                               const int t_end,
                               const int couplstep)
/* print some details about the cvode thermodynamics solver setup
and return a random number generator. Call funciton to set superdroplet
radii to equilibrium wet radius if wetradiiinit is true. */
{
  cvode.print_init_ODEdata(step2dimlesstime(couplstep),
                           step2dimlesstime(t_end));

  for (long unsigned int ii = 0; ii < gridboxes.size(); ++ii)
  {
    set_thermostate(ii, cvode, gridboxes[ii].state);
  }

  if (wetradiiinit)
  {
    set_superdroplets_to_wetradius(gridboxes);
  }

  return std::mt19937(std::random_device()());
}

void set_superdroplets_to_wetradius(std::vector<GridBox> &gridboxes)
/* for each gridbox, set the radius of each superdroplet (SD) to whichever is
larger out of their dry radius or equlibrium wet radius
(given the relative humidity (s_ratio) and temperature of the gridbox).
If relh > maxrelh = 0.95, set each SD's radius to their
equilibrium radius at relh = maxrelh = 0.95 */
{

  const double maxrelh = 0.95;

  for (auto &gbx : gridboxes)
  {
    const double temp = gbx.state.temp;
    const double psat = saturation_pressure(temp);
    const double s_ratio = std::min(maxrelh, supersaturation_ratio(gbx.state.press,
                                                                   gbx.state.qvap, psat));
    for (auto &SDinGBx : gbx.span4SDsinGBx)
    {
      const double equilwetradius = SDinGBx.superdrop.superdroplet_wet_radius(s_ratio, temp);
      const double dryradius = SDinGBx.superdrop.get_dry_radius();
      SDinGBx.superdrop.radius = std::max(dryradius, equilwetradius);
    }
  }
}