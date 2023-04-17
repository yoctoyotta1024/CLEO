// Author: Clara Bayley
// File: timestep_cvodesdm.hpp
/* Header file for functions specifically
to run timstep algoritms for SDM coupled
to Sundials CVODE ODE solver for the
thermodynamics. Coupling can be one-way
or both ways */

#ifndef TIMESTEP_CVODESDM_HPP
#define TIMESTEP_CVODESDM_HPP

#include <random>
#include <vector>

/* Superdroplet Model (SDM) */
#include "sdmgridboxes/run_sdmstep.hpp"
#include "sdmgridboxes/gridbox.hpp"
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/superdropwithgbxindex.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/sdmotion.hpp"
#include "observers/observers.hpp"

/* CVODE thermodynamics ODE solver */
#include "./cvodethermosolver.hpp"

std::vector<ThermoState> start_cvodesdmstep(const Observer auto &observer,
                                            const CvodeThermoSolver &cvode,
                                            std::vector<GridBox> &gridboxes);

int proceedtonext_cvodesdmstep(int t_mdl, const int couplstep,
                               const bool doCouple,
                               const std::vector<ThermoState> &previousstates,
                               const std::vector<GridBox> &gridboxes,
                               CvodeThermoSolver &cvode);

std::vector<ThermoState>
recieve_thermodynamics_from_cvode(const CvodeThermoSolver &cvode,
                                 std::vector<GridBox> &gridboxes);

void send_thermodynamics_to_cvode(const std::vector<ThermoState> &previousstates,
                                      const std::vector<GridBox> &gridboxes,
                                      CvodeThermoSolver &cvode);

inline void set_thermostate(const long unsigned int ii,
                            const CvodeThermoSolver &cvode,
                            ThermoState &state)
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

template <SdMotion M>
void timestep_cvodesdm(const Maps4GridBoxes &gbxmaps,
                       const MoveSuperdropsInDomain<M> &sdmmotion,
                       const SdmProcess auto &sdmprocess,
                       const Observer auto &observer,
                       const bool doCouple,
                       const int t_end,
                       const int couplstep,
                       CvodeThermoSolver &cvode,
                       std::mt19937 &gen,
                       std::vector<GridBox> &gridboxes,
                       std::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* timestep coupled model from t=0 to t=tend. Each coupled step is
length 'outstep' and decomposed into 4 parts: 1) start of step (coupled)
2) run SDM step (independent) 3) run CVODE step (independent)
4) proceed to next step (coupled) */
{
  int t_mdl = 0; // model time is incremented by proceed_tonext_coupledstep
  RunSDMStep<M> sdm();

  while (t_mdl <= t_end)
  {
    /* begin coupled step */
    const std::vector<ThermoState>
        previousstates = start_cvodesdmstep(observer, cvode, gridboxes);

    /* advance SDM by outstep (parallel to CVODE section) */
    sdm.run_sdmstep(t_mdl, couplstep, gbxmaps, sdmmotion,
                sdmprocess, gen, gridboxes, SDsInGBxs);

    /* advance CVODE solver by outstep (parallel to SDM) */
    cvode.run_cvodestep(step2dimlesstime(t_mdl + couplstep));

    /* prepare for next coupled step */
    t_mdl = proceedtonext_cvodesdmstep(t_mdl, couplstep, doCouple,
                                       previousstates, gridboxes, cvode);
  }
}

std::vector<ThermoState> start_cvodesdmstep(const Observer auto &observer,
                                            const CvodeThermoSolver &cvode,
                                            std::vector<GridBox> &gridboxes)
/* communication of thermodynamic state from CVODE solver to SDM.
Sets current thermodynamic state of SDM to match that communicated by
CVODE solver. Then observes each gridbox and then returns vector
of current thermodynamic states (for later use in SDM) */
{
  std::vector<ThermoState>
      currentstates(recieve_thermodynamics_from_cvode(cvode, gridboxes));

  observer.observe_state(gridboxes);

  return currentstates;
}

#endif // TIMESTEP_CVODESDM_HPP