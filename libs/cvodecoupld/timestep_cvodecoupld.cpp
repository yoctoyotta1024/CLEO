// Author: Clara Bayley
// File: timestep_cvodecoupld.hpp
/* Implementation of (non-templated)
functions involved in the timstepping
algoritms for SDM coupled to Sundials
CVODE ODE solver for the thermodynamics.
Coupling is both ways (send and receive) */

#include "./timestep_cvodecoupld.hpp"

std::vector<ThermoState>
  recieve_thermodynamics_from_cvode(const CvodeThermoSolver &cvode,
                                    Kokkos::View<GridBox*> h_gridboxes)
/* get thermo variables from thermodynamics solver and use
these to set ThermoState of each gridbox. Return vector
containing all those Thermostates */
{
  std::vector<ThermoState> currentstates;
  const size_t Ngrid = h_gridboxes.size();
  for (size_t ii(0); ii<Ngrid; ++ii)
  {
    set_thermostate(ii, cvode, h_gridboxes(ii).state);
    currentstates.push_back(h_gridboxes(ii).state);
  }

  return currentstates;
}

int proceedtonext_coupldstep(int t_mdl, const int couplstep,
                             const std::vector<ThermoState> &previousstates,
                             const Kokkos::View<GridBox *> h_gridboxes,
                             CvodeThermoSolver &cvode)
/* sends changes in thermodynamics due to SDM microphysics
to thermodynamics solver (eg. raise in temperature of a
gridbox due to latent heat release).
Then increments timestep by couplstep */
{
  send_thermodynamics_to_cvode(previousstates, h_gridboxes, cvode);

  return t_mdl += couplstep;
}

void send_thermodynamics_to_cvode(const std::vector<ThermoState> &previousstates,
                                  const Kokkos::View<GridBox *> h_gridboxes,
                                  CvodeThermoSolver &cvode)
/* calculate changes in thermodynamics (temp, qv and qc) due to SDM process
affecting ThermoState, then reinitialise cvode solver with those changes */
{
  constexpr int NVARS = 4;
  const size_t Ngrid = h_gridboxes.size();

  std::vector<double> delta_y(Ngrid * NVARS, 0.0);

  for (size_t ii(0); ii<Ngrid; ++ii)
  {
    ThermoState delta_state = gridboxes(ii).state - previousstates.at(ii);

    delta_y[NVARS * ii + 1] = delta_state.temp;
    delta_y[NVARS * ii + 2] = delta_state.qvap;
    delta_y[NVARS * ii + 3] = delta_state.qcond;
  }

  std::vector<double> nodelta(Ngrid * NVARS, 0.0);
  if (delta_y != nodelta)
  {
    cvode.reinitialise(cvode.get_time(), delta_y);
  }
}