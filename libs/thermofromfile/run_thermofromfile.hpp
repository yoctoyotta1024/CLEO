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

void recieve_thermodynamics(const size_t ngbxs, const double time,
                            const ThermodynamicsFromFile &thermodyn,
                            Kokkos::View<GridBox*> h_gridboxes);
/* Sets current thermodynamic state of SDM to match that given
by the ThermodnamicsFromFile 'thermodyn' */

inline Kokkos::Random_XorShift64_Pool<> preparetotimestep()
/* return pool of Kokkos' random number
generators used in SDM */
{
  return Kokkos::Random_XorShift64_Pool<>(std::random_device{}());
}

inline void start_step(const int t_mdl, const int couplstep,
                       const size_t ngbxs,
                       const Observer auto &observer,
                       const ThermodynamicsFromFile &thermodyn,
                       Kokkos::View<GridBox *> h_gridboxes)
/* communication of thermodynamic state
to SDM and observation of SDM gridboxes */
{
  if (t_mdl % couplstep == 0)
  {
    const double time = step2dimlesstime(t_mdl);
    recieve_thermodynamics(ngbxs, time, thermodyn, h_gridboxes);
  }

  if (observer.on_step(t_mdl))
  {
    observer.observe_gridboxes(ngbxs, h_gridboxes);
  }
}

inline int proceedto_next_step(const int t_mdl,
                               const int couplstep,
                               const int obsstep)
/* increments timestep of model to next step. This
is either obsstep or couplstep. Function is also a
placeholder for when communication from SDM to
thermodynamics solver (about changes to thermostates)
could take place if t_mdl was on couplstep. */
{
  const auto nextt = [t_mdl](const int interval)
  {
    return ((t_mdl / interval) + 1) * interval;
  };

  return std::min(nextt(couplstep), nextt(obsstep));
}

template <class MSDs, SdmProcess P, Observer O>
void run_thermofromfile(const Config &config,
                    const RunSDMStep<MSDs, P, O> &sdm,
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
  Kokkos::vector<SuperdropWithGbxindex>
      SDsInGBxs = create_superdrops_from_initSDsfile(config.initSDs_filename,
                                              config.nSDsvec,
                                              config.SDnspace, solute);
  SDsInGBxs.on_device();

  /* vector containing all gridboxes in SDM domain */
  Kokkos::vector<GridBox> gridboxes = create_gridboxes(sdm.gbxmaps, SDsInGBxs);
  gridboxes.on_device(); 
  
  /* prepare model for timestepping */
  auto genpool = preparetotimestep();
  
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
  int t_mdl = 0; // model time is incremented by proceedto_next_step
  
  while (t_mdl <= t_end)
  {
    /* start step (in general involves coupling) */
    gridboxes.on_host(); SDsInGBxs.on_host();
    start_step(t_mdl, couplstep, ngbxs, sdm.observer,
               thermodyn, gridboxes.view_host());

    /* advance SDM by couplstep
    (optionally concurrent to thermodynamics solver) */
    gridboxes.on_device(); SDsInGBxs.on_device();
    sdm.run_sdmstep(t_mdl, couplstep, genpool, gridboxes, SDsInGBxs);

    /* advance thermodynamics solver
    (optionally concurrent to SDM) */
    thermodyn.run_thermostep(t_mdl, couplstep);

    /* proceed to next step (in general involves coupling) */
    t_mdl = proceedto_next_step(t_mdl, couplstep,
                                sdm.observer.get_interval());
  }
}

#endif // RUN_THERMOFROMFILE_HPP