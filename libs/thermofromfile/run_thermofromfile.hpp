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
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/sdmotion.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "observers/observers.hpp"

/* Thermodynamics Solver */
#include "thermofromfile/thermodynamicsfromfile.hpp"

namespace dlc = dimless_constants;

std::mt19937 preparetotimestep()
{
  std::cout << "any preparation? \n";

  return std::mt19937(std::random_device()());
}

void start_step()
{
  std::cout << "start step\n";
}

int proceedto_next_step(int t_mdl, const int couplstep)
{
  std::cout << "at the very least, advance to next step\n";
  return t_mdl + couplstep;
}

void timestep_thermofromfile(const int t_end,
                            const int couplstep,
                            const RunSDMStep<auto, auto, auto> &sdm,
                            const ThermodynamicsFromFile &thermodyn,
                            std::mt19937 &gen,
                            std::vector<GridBox> &gridboxes,
                            std::vector<SuperdropWithGbxindex> &SDsInGBxs);

void run_thermofromfile(const Config &config,
                    const RunSDMStep<auto, auto, auto> &sdm,
                    const int t_end, const int couplstep)
/* create superdroplets and gridboxes and then run uncoupled 
superdroplet model (SDM) using thermodynamics read from files */
{
  /* create thermodynamics from file */
  const unsigned int ngridboxes = sdm.ngridboxes;
  const ThermodynamicsFromFile thermodyn(ngridboxes);
  
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
    start_step();

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