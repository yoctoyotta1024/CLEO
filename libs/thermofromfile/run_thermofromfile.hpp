// Author: Clara Bayley
// File: run_thermofromfile.hpp
/* Header file for functions specifically
to run uncoupled SDM where thermodynamics
are read from file */

#ifndef RUN_THERMOFROMFILE_HPP
#define RUN_THERMOFROMFILE_HPP

#include <random>
#include <vector>
#include <iostream>
#include <memory>

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
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/sdmotion.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "observers/observers.hpp"

/* Thermodynamics Solver */
#include "./thermodynamicsfromfile.hpp"

namespace dlc = dimless_constants;

void print_state(const std::vector<GridBox> &gridboxes)
{
  for (long unsigned int ii = 0; ii < gridboxes.size(); ++ii)
  {
    std::cout << "gbx " << ii << ", " << gridboxes[ii].state.press << ", ";
    std::cout << gridboxes[ii].state.temp<< ", ";
    std::cout << gridboxes[ii].state.qvap << ", ";
    std::cout << gridboxes[ii].state.qcond << ", ";
    std::cout << gridboxes[ii].state.wvel << "\n";
  }
}

void timestep_thermofromfile(const int t_end,
                            const int couplstep,
                            const RunSDMStep<auto, auto, auto> &sdm,
                            const ThermodynamicsFromFile &thermodyn,
                            std::mt19937 &gen,
                            std::vector<GridBox> &gridboxes,
                            std::vector<SuperdropWithGbxindex> &SDsInGBxs);

void recieve_thermodynamics(const int time,
                            const ThermodynamicsFromFile &thermodyn,
                            std::vector<GridBox> &gridboxes);
/* Sets current thermodynamic state of SDM to match that given
by the ThermodnamicsFromFile 'thermodyn' */

inline std::mt19937 preparetotimestep()
/* return random number generator used in SDM */
{
  return std::mt19937(std::random_device()());
}

inline void start_step(const int t_mdl,
                       const Observer auto &observer,
                       const ThermodynamicsFromFile &thermodyn,
                       std::vector<GridBox> &gridboxes)
/* communication of thermodynamic state
to SDM and observation of SDM gridboxes */
{
  recieve_thermodynamics(t_mdl, thermodyn, gridboxes);

  observer.observe_state(gridboxes);
}

inline int proceedto_next_step(int t_mdl, const int couplstep)
/* increments timestep of model by couplstep.
This function is also placeholder for moment when 
communication from SDM to thermodynamic solver
about thermodynamic state (changes) is possible. */
{
  return t_mdl + couplstep;
}

void run_thermofromfile(const Config &config,
                    const RunSDMStep<auto, auto, auto> &sdm,
                    const int t_end, const int couplstep)
/* create superdroplets and gridboxes and then run uncoupled 
superdroplet model (SDM) using thermodynamics read from files */
{
  /* create thermodynamics from file */
  const ThermodynamicsFromFile thermodyn(config);
  
  /* vector containing all superdroplets within a
  struct that also holds their associated gridbox index.
  (all superdroplets have same solute properties) */
  const auto solute(std::make_shared<const SoluteProperties>());
  std::vector<SuperdropWithGbxindex>
      SDsInGBxs = create_superdrops_from_initSDsfile(config.initSDs_filename,
                                              config.nSDsvec,
                                              config.SDnspace, solute);

  /* vector containing all gridboxes in SDM domain */
  std::vector<GridBox> gridboxes = create_gridboxes(sdm.gbxmaps, SDsInGBxs);

  /* prepare model for timestepping */
  auto gen = preparetotimestep();

  /* run model from t=0 to t=t_end */
  timestep_thermofromfile(t_end, couplstep, sdm, thermodyn,
                          gen, gridboxes, SDsInGBxs);

  std::cout << "\n ---- Uncoupled SDM Run Complete ---- \n";
}

void timestep_thermofromfile(const int t_end,
                            const int couplstep,
                            const RunSDMStep<auto, auto, auto> &sdm,
                            const ThermodynamicsFromFile &thermodyn,
                            std::mt19937 &gen,
                            std::vector<GridBox> &gridboxes,
                            std::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* timestep model from t=0 to t=tend. Each step is
length 'couplstep' and is decomposed into 4 parts:
1) start of step (optionally coupled)
2) run SDM step (independent, optionally concurrent)
3) run thermodynamics (independent, optionally concurrent)
4) proceed to next step (optionally coupled) */
{
  int t_mdl = 0; // model time is incremented by proceedto_next_step
  
  while (t_mdl <= t_end)
  {
    /* start step (in general involves coupling) */
    start_step(t_mdl, sdm.observer, thermodyn, gridboxes);

    /* advance SDM by couplstep
    (optionally concurrent to thermodynamics solver) */
    sdm.run_sdmstep(t_mdl, couplstep, gen, gridboxes, SDsInGBxs);

    /* advance thermodynamics solver by couplstep
    (optionally concurrent to SDM) */
    thermodyn.run_thermostep(couplstep);

    /* proceed to next step (in general involves coupling) */
    t_mdl = proceedto_next_step(t_mdl, couplstep);
  }
}

#endif // RUN_THERMOFROMFILE_HPP