// Author: Clara Bayley
// File: run_thermofromfile.hpp
/* Header file for functions specifically
to run uncoupled SDM where thermodynamics
are read from file */

#ifndef RUN_THERMOFROMFILE_HPP
#define RUN_THERMOFROMFILE_HPP

#include <random>
#include <iostream>
#include <memory>
#include <cmath>

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>
#include <Kokkos_Random.hpp>

/* Coupled model setup */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

/* Superdroplet Model (SDM) */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/superdropwithgbxindex.hpp"
#include "sdmgridboxes/movesuperdropsindomain.hpp"
#include "sdmgridboxes/gridbox.hpp"
#include "sdmgridboxes/sdmtimesteps.hpp"
#include "sdmgridboxes/runsdmstep.hpp"
#include "sdmgridboxes/sdmotion.hpp"
#include "sdmgridboxes/detectors.hpp"
#include "sdmgridboxes/detectors_ptr.hpp"
#include "sdmgridboxes/logbooks.hpp"
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "observers/observers.hpp"

/* Thermodynamics Solver */
#include "./thermodynamicsfromfile.hpp"

namespace dlc = dimless_constants;

template <class MSDs, SdmProcess P, Observer O>
void timestep_thermofromfile(const int t_end,
                            const int couplstep,
                            const RunSDMStep<MSDs, P, O> &sdm,
                            ThermodynamicsFromFile &thermodyn,
                            Kokkos::Random_XorShift64_Pool<> &genpool,
                            Kokkos::vector<GridBox> &gridboxes,
                            Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs);

void receive_thermodynamics(const int t_mdl, const int couplstep,
                            const size_t ngbxs,
                            const ThermodynamicsFromFile &thermodyn,
                            Kokkos::View<GridBox *> h_gridboxes);
/* updates time in each gbx thermodynamic state to
match t_mdl and receives thermodynamics from thermodyanmic
solver 'thermodyn' if on couplstep */

void receive_thermodynamics_from_thermodyn(const size_t ngbxs,
                            const ThermodynamicsFromFile &thermodyn,
                            Kokkos::View<GridBox*> h_gridboxes);
/* Sets current thermodynamic state of SDM to match that given
by the ThermodnamicsFromFile 'thermodyn' */

template <class MSDs, SdmProcess P, Observer O>
inline Kokkos::Random_XorShift64_Pool<>
preparetotimestep(const RunSDMStep<MSDs, P, O> &sdm)
/* return pool of Kokkos' random number
generators used in SDM and prepare observer */
{
  sdm.observer.prepare(sdm.logbooks);

  return Kokkos::Random_XorShift64_Pool<>(std::random_device{}());
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

inline int start_step(const int t_mdl, const int couplstep,
                       const size_t ngbxs,
                       const Observer auto &observer,
                       const DetectorLogbooks &logbooks,
                       const ThermodynamicsFromFile &thermodyn,
                       Kokkos::View<GridBox *> h_gridboxes)
/* optional communication of thermodynamic state
to SDM and observation of SDM gridboxes. returns step size 
to take given current time t_mdl */
{
  receive_thermodynamics(t_mdl, couplstep, ngbxs,
                         thermodyn, h_gridboxes);

  if (observer.on_step(t_mdl))
  {
    observer.observe(ngbxs, h_gridboxes, logbooks);
  }

  return stepsize(t_mdl, couplstep, observer.get_interval());
}

inline int proceedto_next_step(const int t_mdl, const int onestep)
/* returns incremented timestep 't_mdl' of model
by 'onestep'. Function is also a placeholder for
when communication from SDM to thermodynamics
solver (about changes to thermostates) could
take place if t_mdl was on couplstep. */
{
  return t_mdl + onestep;
}

template <class MSDs, SdmProcess P, Observer O>
void run_thermofromfile(const Config &config,
                        const RunSDMStep<MSDs, P, O> &sdm,
                        const CreateDetectorsPtr auto &dtrs,
                        const int t_end, const int couplstep)
/* create superdroplets and gridboxes and then run uncoupled
superdroplet model (SDM) using thermodynamics read from files */
{
  Kokkos::Timer kokkostimer;

  /* create thermodynamics from file */
  const size_t nsteps = ceil(t_end / couplstep) + 1;
  ThermodynamicsFromFile thermodyn(config, sdm.gbxmaps.ndims, nsteps);

  /* vector containing all superdroplets within a
  struct that also holds their associated gridbox index.
  (all superdroplets have same solute properties) */
  const auto solute(std::make_shared<const SoluteProperties>());
  Kokkos::vector<SuperdropWithGbxindex> SDsInGBxs =
      create_superdrops_from_initSDsfile(config.initSDs_filename,
                                         config.nSDsvec,
                                         config.SDnspace, solute);
  SDsInGBxs.on_device();

  /* vector containing all gridboxes in SDM domain */
  Kokkos::vector<GridBox> gridboxes =
      create_gridboxes(sdm.gbxmaps, sdm.logbooks, dtrs, SDsInGBxs);
  gridboxes.on_device(); 
  
  /* prepare model for timestepping */
  auto genpool = preparetotimestep(sdm);
  
  const double t1 = kokkostimer.seconds();
  /* run model from t=0 to t=t_end */
  timestep_thermofromfile(t_end, couplstep, sdm, thermodyn,
                          genpool, gridboxes, SDsInGBxs);
  const double t2 = kokkostimer.seconds();
  
  std::cout << "\n ---- Uncoupled SDM Run Complete ---- \n"
            << "       Duration: " << t2 << "s ----- \n"
            << "       Initialisation: " << t1 << "s ----- \n"
            << "       Timestepping: " << t2 - t1 << "s ----- \n"
            << "\n ------------------------------------ \n";
}

template <class MSDs, SdmProcess P, Observer O>
void timestep_thermofromfile(const int t_end,
                            const int couplstep,
                            const RunSDMStep<MSDs, P, O> &sdm,
                            ThermodynamicsFromFile &thermodyn,
                            Kokkos::Random_XorShift64_Pool<> &genpool,
                            Kokkos::vector<GridBox> &gridboxes,
                            Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* timestep model from t=0 to t=tend. Each step is
length 'couplstep' and is decomposed into 4 parts:
1) start of step (optionally coupled)
2) run SDM step (independent, optionally concurrent)
3) run thermodynamics (independent, optionally concurrent)
4) proceed to next step (optionally coupled) */
{
  const size_t ngbxs(gridboxes.size());
  int t_mdl = 0; // model time is incremented at proceedto_next_step
  
  while (t_mdl <= t_end)
  {
    /* start step (in general involves coupling) */
    gridboxes.on_host(); SDsInGBxs.on_host();
    const int onestep = start_step(t_mdl, couplstep, ngbxs,
                                   sdm.observer, sdm.logbooks,
                                   thermodyn, gridboxes.view_host());

    /* advance SDM from t_mdl to t_mdl + onestep
    (optionally concurrent to thermodynamics solver) */
    gridboxes.on_device(); SDsInGBxs.on_device();
    sdm.run_sdmstep(t_mdl, onestep, genpool, gridboxes, SDsInGBxs);

    /* advance thermodynamics solver
    (optionally concurrent to SDM) */
    thermodyn.run_thermostep(t_mdl, couplstep);

    /* proceed to next step (in general involves coupling) */
    t_mdl = proceedto_next_step(t_mdl, onestep);
  }
}

#endif // RUN_THERMOFROMFILE_HPP