// Author: Clara Bayley
// File: run_cvodecoupld.cpp
/* Functions involved specifically
in running SDM coupled to Sundials
CVODE ODE solver for the thermodynamics.
Coupling is both ways (send and receive) */

#include "./run_cvodecoupld.hpp"

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

Kokkos::Random_XorShift64_Pool<>
preparetotimestep(CvodeThermoSolver &cvode,
                  std::vector<GridBox> &gridboxes,
                  const bool wetradiiinit,
                  const int t_end,
                  const int couplstep)
/* print some details about the cvode thermodynamics solver setup
and return pool of Kokkos' random number generator. Call
function to set superdroplet radii to equilibrium wet radius
if wetradiiinit is true. */
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

  return Kokkos::Random_XorShift64_Pool<>(std::random_device{}()); // pool of Kokkos' random number generators
}