// Author: Clara Bayley
// File: timestep_cvodecoupld.hpp
/* Implementation of (non-templated)
functions involved in the timstepping
algoritms for SDM coupled to Sundials
CVODE ODE solver for the thermodynamics.
Coupling can be one-way or both ways */

#include "./timestep_cvodecoupld.hpp"

std::vector<ThermoState>
  recieve_thermodynamics_from_cvode(const CvodeThermoSolver &cvode,
                                 std::vector<GridBox> &gridboxes)
/* get thermo variables from thermodynamics solver and use
these to set ThermoState of each gridbox. Return vector
containing all those Thermostates */
{
  std::vector<ThermoState> currentstates;
  for (long unsigned int ii = 0; ii < gridboxes.size(); ++ii)
  {
    set_thermostate(ii, cvode, gridboxes[ii].state);
    currentstates.push_back(gridboxes[ii].state);
  }

  return currentstates;
}

int proceedtonext_coupldstep(int t_mdl, const int couplstep,
                               const std::vector<ThermoState> &previousstates,
                               const std::vector<GridBox> &gridboxes,
                               CvodeThermoSolver &cvode)
/* sends changes in thermodynamics due to SDM microphysics
to thermodynamics solver (eg. raise in temperature of a
gridbox due to latent heat release) */
{

  send_thermodynamics_to_cvode(previousstates, gridboxes, cvode);

  return t_mdl += couplstep;
}

void send_thermodynamics_to_cvode(const std::vector<ThermoState> &previousstates,
                                      const std::vector<GridBox> &gridboxes,
                                      CvodeThermoSolver &cvode)
/* calculate changes in thermodynamics (temp, qv and qc) due to SDM process
affecting ThermoState, then reinitialise cvode solver with those changes */
{
  constexpr int NVARS = 4;
  const long unsigned int gridN = gridboxes.size();

  std::vector<double> delta_y(gridN * NVARS, 0.0);

  for (long unsigned int ii = 0; ii < gridN; ++ii)
  {
    ThermoState delta_state = gridboxes[ii].state - previousstates[ii];

    delta_y[NVARS * ii + 1] = delta_state.temp;
    delta_y[NVARS * ii + 2] = delta_state.qvap;
    delta_y[NVARS * ii + 3] = delta_state.qcond;
  }

  std::vector<double> nodelta(gridN * NVARS, 0.0);
  if (delta_y != nodelta)
  {
    cvode.reinitialise(cvode.get_time(), delta_y);
  }
}