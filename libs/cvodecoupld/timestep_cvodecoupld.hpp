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

std::vector<ThermoState> start_coupldstep(const size_t ngbxs,
                                          const Observer auto &observer,
                                          const CvodeThermoSolver &cvode,
                                          const Kokkos::View<GridBox *> h_gridboxes);

int proceedtonext_coupldstep(int t_mdl, const int couplstep,
                             const size_t ngbxs,
                             const std::vector<ThermoState> &previousstates,
                             const Kokkos::View<GridBox *> h_gridboxes,
                             CvodeThermoSolver &cvode);

std::vector<ThermoState>
recieve_thermodynamics_from_cvode(const size_t ngbxs,
                                  const CvodeThermoSolver &cvode,
                                  Kokkos::View<GridBox *> h_gridboxes);

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
  int t_mdl = 0; // model time is incremented by proceedtonext_coupledstep
  
  while (t_mdl <= t_end)
  {
    /* start step (in general involves coupling) */
    const std::vector<ThermoState>
        previousstates = start_coupldstep(ngbxs,
                                          sdm.observer, cvode,
                                          gridboxes.view_host());

    /* advance SDM by couplstep (optionally
    concurrent to CVODE thermodynamics solver) */
    gridboxes.on_device(); SDsInGBxs.on_device();
    sdm.run_sdmstep(t_mdl, couplstep, genpool, gridboxes, SDsInGBxs);

    /* advance CVODE thermodynamics solver by
    couplstep (optionally concurrent to SDM) */
    cvode.run_cvodestep(step2dimlesstime(t_mdl + couplstep));

    /* prepare for next coupled step (in general involves coupling) */
    gridboxes.on_host(); SDsInGBxs.on_host();
    t_mdl = proceedtonext_coupldstep(t_mdl, couplstep, ngbxs,
                                     previousstates,
                                     gridboxes.view_host(),
                                     cvode);
  }
}

std::vector<ThermoState> start_coupldstep(const size_t ngbxs,
                                          const Observer auto &observer,
                                          const CvodeThermoSolver &cvode,
                                          Kokkos::View<GridBox *> h_gridboxes)
/* communication of thermodynamic state from CVODE solver to SDM.
Sets current thermodynamic state of SDM to match that communicated by
CVODE solver. Then observes each gridbox and then returns vector
of current thermodynamic states (for later use in SDM) */
{
  std::vector<ThermoState>
      currentstates(recieve_thermodynamics_from_cvode(ngbxs, cvode, h_gridboxes));

  observer.observe_gridboxes(ngbxs, h_gridboxes);

  return currentstates;
}

#endif // TIMESTEP_CVODECOUPLD_HPP