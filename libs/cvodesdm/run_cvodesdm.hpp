// Author: Clara Bayley
// File: run_cvodesdm.hpp
/* Header file for functions specifically
to run SDM coupled to Sundials CVODE ODE
solver for the thermodynamics.
Coupling can be one-way or both ways */

#ifndef RUN_CVODESDM_HPP
#define RUN_CVODESDM_HPP

#include <random>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>

/* Coupled model setup */
#include "claras_SDconstants.hpp"
#include "./timestep_cvodesdm.hpp"
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

/* CVODE thermodynamics ODE solver */
#include "./cvodethermosolver.hpp"

namespace dlc = dimless_constants;

std::mt19937 preparetotimestep(CvodeThermoSolver &cvode,
                                  std::vector<GridBox> &gridboxes,
                                  const bool wetradiiinit,
                                  const int t_end,
                                  const int couplstep);

std::vector<double> initcvodethermo(const size_t num_gridboxes,
                              const Config &config);

void set_superdroplets_to_wetradius(std::vector<GridBox> &gridboxes);

void run_cvodesdm(const Config &config,
                  const RunSDMStep<auto, auto, auto> &sdm,
                  const int t_end, const int couplstep)
/* create CVODE thermodynamics solver, superdroplets and gridboxes and
then run superdroplet model (SDM) coupled to the thermodynamics solver */
{
  /* CVODE thermodynamics solver */
  const unsigned int ngridboxes = sdm.ngridboxes;
  CvodeThermoSolver cvode(config, initcvodethermo(ngridboxes, config));

  /* vector containing all superdroplets within a struct that also holds their
  associated gridbox index. (all superdroplets have same solute properties) */
  const auto solute(std::make_shared<const SoluteProperties>());
  std::vector<SuperdropWithGbxindex>
      SDsInGBxs = create_superdrops_from_initSDsfile(config.initSDs_filename,
                                              config.nSDsvec,
                                              config.SDnspace, solute);

  /* vector containing all gridboxes that makeup the SDM domain */
  std::vector<GridBox> gridboxes = create_gridboxes(sdm.gbxmaps, SDsInGBxs);

  /* prepare coupled model for timestepping */
  auto gen = preparetotimestep(cvode, gridboxes, config.wetradiiinit,
                                  t_end, couplstep);

  /* run coupled model from t=0 to t=t_end */
  timestep_cvodesdm(t_end, couplstep, config.doAlterThermo, sdm, cvode,
                    gen, gridboxes, SDsInGBxs);

  std::cout << "\n ---- CVODE-SDM Coupled Model Complete ---- \n";
}

#endif // RUN_CVODESDM_HPP