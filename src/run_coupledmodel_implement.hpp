// Author: Clara Bayley
// File: run_coupledmodel.hpp
/* continuation of 'run_coupledmodel.hpp',
included at start of that file. this seperate
file just contains start of that file but is
seperated to increase readability */

#ifndef RUN_COUPLEDMODEL_IMPLEMENT_HPP
#define RUN_COUPLEDMODEL_IMPLEMENT_HPP

#include <random>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>

/* constants and coupled model setup */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"
#include "./timesteps.hpp"

/* sdm gridboxes setup */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/superdrops_in_gridboxes.hpp"
#include "sdmgridboxes/movement_in_domain.hpp"
#include "sdmgridboxes/gridbox.hpp"
#include "observers/observers.hpp"

/* Superdroplet Model (SDM) files */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/sdmmotion.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/thermostate.hpp"

/* CVODE ODE thermodynamics solver files */
#include "thermo_solver/cvodethermosolver.hpp"

namespace dlc = dimless_constants;

/* ----------- implementation in run_coupledmodel.hpp ----------- */
void run_cvodeSDM_coupledmodel(const Config &config,
                               const Timesteps &mdlsteps,
                               const Maps4GridBoxes &mdlmaps,
                               const SdmProcess auto &sdmprocess,
                               const SdmMotion &sdmmotion,
                               const Observer auto &observer);
/* create CVODE thermodynamics solver, superdroplets and gridboxes and
then run superdroplet model (SDM) coupled to the thermodynamics solver */

void timestep_coupledmodel(const Timesteps &mdlsteps,
                           const Maps4GridBoxes &mdlmaps,
                           const SdmProcess auto &sdmprocess,
                           const SdmMotion &sdmmotion,
                           const Observer auto &observer,
                           const bool doCouple,
                           CvodeThermoSolver &cvode,
                           std::mt19937 &gen,
                           std::vector<GridBox> &gridboxes,
                           std::vector<SuperdropWithGbxindex> &SDsInGBxs);
/* timestep coupled model from t=0 to t=tend. Each coupled step is
length 'outstep' and decomposed into 4 parts: 1) start of step (coupled)
2) run SDM step (independent) 3) run CVODE step (independent)
4) proceed to next step (coupled) */
/* -------------------------------------------------------------- */

/* ----------- implementation in run_coupledmodel.cpp ----------- */
std::vector<double> init_thermodynamics(const size_t num_gridboxes,
                                        const Config &config);
/* return vector of dimensionless initial conditions
for thermodyanmic variables (p, temp, qv, qc) to
initialise cvode thermodynamics solver */

std::mt19937 prepare_coupledmodel(const Timesteps &mdlsteps, CvodeThermoSolver &cvode,
                                  std::vector<GridBox> &gridboxes,
                                  const bool wetradiiinit);
/* print some details about the cvode thermodynamics solver setup and
return a random number generator. Call funciton to set superdroplet radii
to equilibrium wet radius if wetradiiinit is true. */

std::vector<ThermoState> set_thermodynamics_from_cvodesolver(std::vector<GridBox> &gridboxes,
                                                             const CvodeThermoSolver &cvode);
/* get thermo variables from thermodynamics solver and use
these to set ThermoState of each gridbox. Return vector
containing all those Thermostates */

int proceed_tonext_coupledstep(int t_out, const int outstep,
                               const bool doCouple,
                               const std::vector<ThermoState> &previousstates,
                               std::vector<GridBox> &gridboxes,
                               CvodeThermoSolver &cvode);
/* exchanges superdroplets between gridboxes and sends
changes in thermodynamics due to SDM microphysics to thermodynamics solver
(eg. raise in temperature of a gridbox due to latent heat release) */
/* -------------------------------------------------------------- */

inline void printfinish_coupledmodel()
/* print statement declaring coupled model completed*/
{
  std::cout << "\n ---- Coupled Model Run Complete ---- \n";
}

inline void set_thermostate(const long unsigned int ii,
                            ThermoState &state,
                            const CvodeThermoSolver &cvode)
/* set values of the ThermoState instance's members (time,
p, temp, qv, qc, etc.) using data sent from the
thermodyanics ODE solver (cvode) */
{
  state.time = cvode.get_time();
  state.press = cvode.get_pressure(ii);
  state.temp = cvode.get_temperature(ii);
  state.qvap = cvode.get_qvap(ii);
  state.qcond = cvode.get_qcond(ii);
}

inline int exchange_or_outstep(const int t_out, const int outstep,
                               const int xchangestep)
/* given current time, t_out, work out which event (exchange or output)
is next to occur and return the time of the sooner event */
{
  const int next_xchange = ((t_out / xchangestep) + 1) * xchangestep; // t of next xchange
  const int next_out = ((t_out / outstep) + 1) * outstep;             // t of next output

  return std::min(next_xchange, next_out);
}

inline void exchanges_between_gridboxes(const Maps4GridBoxes &mdlmaps,
                                        const SdmMotion &sdmmotion,
                                        std::vector<SuperdropWithGbxindex> &SDsInGBxs,
                                        std::vector<GridBox> &gridboxes)
{
  move_superdrops_in_domain(mdlmaps, sdmmotion, SDsInGBxs, gridboxes);
}

std::vector<ThermoState> start_coupledstep(const Observer auto &observer,
                                           std::vector<GridBox> &gridboxes,
                                           const CvodeThermoSolver &cvode)
/* communication of thermodynamic state from CVODE solver to SDM.
Sets current thermodynamic state of SDM to match that communicated by
CVODE solver. Then observes each gridbox and then returns vector
of current thermodynamic states (for later use in SDM) */
{
  std::vector<ThermoState> currentstates(
      set_thermodynamics_from_cvodesolver(gridboxes, cvode));

  observer.observe_state(gridboxes);

  return currentstates;
}

void run_sdmstep(const int t_out, const int outstep,
                 const int xchangestep,
                 const SdmProcess auto &sdmprocess,
                 const SdmMotion &sdmmotion,
                 const Maps4GridBoxes &mdlmaps,
                 std::mt19937 &gen,
                 std::vector<GridBox> &gridboxes,
                 std::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* run SDM for each gridbox from time t_out to t_out+outstep
with subtimestepping such that each output timestep (outstep)
can be subdivided to allow the exchange of superdroplets between
gridboxes and the SDM process to occur at smaller time intervals */
{

  int t = t_out;

  while (t < t_out + outstep)
  {
    /* nextt is t of next exchange and/or t of next outstep */
    const int nextt = exchange_or_outstep(t, outstep, xchangestep);

    /* run SDM process for all gridboxes from t_out to nextt
    using SDM subt timestepping routine */
    for (auto &gbx : gridboxes)
    {
      for (int subt = t; subt < nextt;
           subt = sdmprocess.next_step(subt))
      {
        sdmprocess.run_step(subt, gbx.span4SDsinGBx,
                            gbx.state, gen);
      }
    }

    /* do exchange if timestep is on exchange event */
    if (nextt % xchangestep == 0)
    {
      exchanges_between_gridboxes(mdlmaps, sdmmotion,
                                  SDsInGBxs, gridboxes);
    }

    t = nextt;
  }
}

#endif // RUN_COUPLEDMODEL_IMPLEMENT_HPP