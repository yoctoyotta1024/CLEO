// Author: Clara Bayley
// File: run_cvodecoupld.cpp
/* Functions involved specifically
in running SDM coupled to Sundials
CVODE ODE solver for the thermodynamics.
Coupling is both ways (send and receive) */

#include "./run_cvodecoupld.hpp"

std::vector<double> initcvodethermo(const size_t ngridboxes,
                                        const Config &config)
/* return vector of dimensionless initial conditions
for thermodyanmic variables (p, temp, qv, qc) to
initialise cvode thermodynamics solver */
{
  constexpr int NVARS = 4;                  // no. (distinct) variables per grid box
  const size_t neq = NVARS * ngridboxes; // total no. variables

  const double p_init = config.P_INIT / dlc::P0;
  const double temp_init = config.TEMP_INIT / dlc::TEMP0;
  const double vapourp_init = saturation_pressure(temp_init) * config.relh_init / 100.0;
  const double qv_init = vapourpressure_2_massmixratio(vapourp_init, p_init);
  const double qc_init = config.qc_init;

  std::vector<double> y_init(neq);
  for (size_t k = 0; k < neq; k += NVARS)
  {
    y_init.at(k) = p_init;
    y_init.at(k + 1) = temp_init;
    y_init.at(k + 2) = qv_init;
    y_init.at(k + 3) = qc_init;
  }

  return y_init;
}