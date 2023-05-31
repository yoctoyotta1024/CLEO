// Author: Clara Bayley
// File: timestep_cvodecoupld.hpp
/* Header file for functions specifically
to run timstep algoritms for SDM coupled
to Sundials CVODE ODE solver for the
thermodynamics. Coupling is in general 
both ways (send and receive) */

#ifndef TIMESTEP_CVODECOUPLD_HPP
#define TIMESTEP_CVODECOUPLD_HPP

#include <random>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Vector.hpp>

/* Superdroplet Model (SDM) */
#include "sdmgridboxes/runsdmstep.hpp"
#include "sdmgridboxes/gridbox.hpp"
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/sdmtimesteps.hpp"
#include "sdmgridboxes/superdropwithgbxindex.hpp"
#include "sdmgridboxes/sdmotion.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "observers/observers.hpp"

/* CVODE thermodynamics ODE solver */
#include "./cvodethermosolver.hpp"

std::vector<ThermoState>
receive_thermodynamics(const int t_mdl,
                       const int couplstep,
                       const size_t ngbxs,
                       const CvodeThermoSolver &cvode,
                       Kokkos::View<GridBox *> h_gridboxes);

std::vector<ThermoState>
recieve_thermodynamics_from_cvode(const size_t ngbxs,
                                  const CvodeThermoSolver &cvode,
                                  Kokkos::View<GridBox *> h_gridboxes);

int proceedto_next_step(const int t_mdl,
                        const int onestep,
                        const int couplstep,
                        const size_t ngbxs,
                        const std::vector<ThermoState> &previousstates,
                        const Kokkos::View<GridBox *> h_gridboxes,
                        CvodeThermoSolver &cvode);

void send_thermodynamics_to_cvode(const size_t ngbxs,
                                  const std::vector<ThermoState> &previousstates,
                                  const Kokkos::View<GridBox *> h_gridboxes,
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

inline int stepsize(const int t_mdl, const int couplstep, const int obsstep)
/* returns size of next step of model ('onestep')
given current time t_mdl, so that next time
(t_next = t_mdl + onestep) is time of obs or coupl */
{
  const auto next_step = [t_mdl](const int interval)
  {
    return ((t_mdl / interval) + 1) * interval;
  };

  /* t_next is smaller out of time of next coupl and obs */
  const int next_coupl(next_step(couplstep));
  const int next_obs(next_step(obsstep));

  return std::min(next_coupl, next_obs) - t_mdl;
}

inline std::vector<ThermoState>
start_coupldstep(const int t_mdl,
                 const int couplstep,
                 const size_t ngbxs,
                 const Observer auto &observer,
                 const CvodeThermoSolver &cvode,
                 Kokkos::View<GridBox *> h_gridboxes)
/* communication of thermodynamic state from CVODE solver to SDM.
Sets current thermodynamic state of SDM to match that communicated by
CVODE solver. Then observes each gridbox and then returns vector
of current thermodynamic states (for later use in SDM) */
{
  auto currentstates = receive_thermodynamics(t_mdl, couplstep,
                                         ngbxs, cvode, h_gridboxes);

  if (observer.on_step(t_mdl))
  {
    observer.observe_gridboxes(ngbxs, h_gridboxes);
  }
  
  return currentstates;
}

void timestep_cvodecoupld(const int t_end,
                          const int couplstep,
                          const RunSDMStep<auto, auto, auto> &sdm,
                          CvodeThermoSolver &cvode,
                          Kokkos::Random_XorShift64_Pool<> &genpool,
                          Kokkos::vector<GridBox> &gridboxes,
                          Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* timestep coupled model from t=0 to t=tend. Each coupled step is
length 'couplstep' and decomposed into 4 parts:
1) start of step (coupled)
2) run SDM step (independent, optionally concurrent)
3) run CVODE step (independent, optionally concurrent)
4) proceed to next step (coupled) */
{
  const size_t ngbxs(gridboxes.size());
  int t_mdl = 0; // model time is incremented at proceedto_next_step
  
  while (t_mdl <= t_end)
  {
    const int onestep = stepsize(t_mdl, couplstep,
                                sdm.observer.get_interval());

    /* start step (in general involves coupling) */
    const std::vector<ThermoState>
        previousstates = start_coupldstep(t_mdl, couplstep, ngbxs,
                                          sdm.observer, cvode,
                                          gridboxes.view_host());

    /* advance SDM by couplstep (optionally
    concurrent to CVODE thermodynamics solver) */
    gridboxes.on_device(); SDsInGBxs.on_device();
    sdm.run_sdmstep(t_mdl, onestep, genpool, gridboxes, SDsInGBxs);

    /* advance CVODE thermodynamics solver to one
    coupled step 'couplstep' (optionally concurrent to SDM) */
    cvode.run_cvodestep(t_mdl, couplstep,
                        step2dimlesstime(t_mdl + couplstep));

    /* prepare for next coupled step (in general involves coupling) */
    gridboxes.on_host(); SDsInGBxs.on_host();
    t_mdl = proceedto_next_step(t_mdl, onestep, couplstep,
                                ngbxs,
                                previousstates,
                                gridboxes.view_host(),
                                cvode);
  }
}

#endif // TIMESTEP_CVODECOUPLD_HPP