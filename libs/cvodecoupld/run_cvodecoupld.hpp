// Author: Clara Bayley
// File: run_cvodecoupld.hpp
/* Header file for functions specifically
to run SDM coupled to Sundials CVODE ODE
solver for the thermodynamics.
Coupling is both ways (send and receive) */

#ifndef RUN_CVODECOUPLD_HPP
#define RUN_CVODECOUPLD_HPP

#include <random>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>
#include <Kokkos_Random.hpp>

/* Coupled model setup */
#include "claras_SDconstants.hpp"
#include "./timestep_cvodecoupld.hpp"
#include "initialisation/config.hpp"

/* Superdroplet Model (SDM) */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/superdropwithgbxindex.hpp"
#include "sdmgridboxes/movesuperdropsindomain.hpp"
#include "sdmgridboxes/gridbox.hpp"
#include "sdmgridboxes/sdmtimesteps.hpp"
#include "sdmgridboxes/sdmotion.hpp"
#include "sdmgridboxes/detectors.hpp"
#include "sdmgridboxes/detectors_ptr.hpp"
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "observers/observers.hpp"

/* CVODE thermodynamics ODE solver */
#include "./cvodethermosolver.hpp"

namespace dlc = dimless_constants;

template <class MSDs, SdmProcess P, Observer O>
Kokkos::Random_XorShift64_Pool<>
preparetotimestep(const RunSDMStep<MSDs, P, O> &sdm,
                  CvodeThermoSolver &cvode,
                  Kokkos::vector<GridBox> &gridboxes,
                  const bool wetradiiinit,
                  const int t_end,
                  const int couplstep);
/* print some details about the cvode thermodynamics solver setup
and return pool of Kokkos' random number generator. Call
function to set superdroplet radii to equilibrium wet radius
if wetradiiinit is true. Call sdm observer's prepare function */

std::vector<double> initcvodethermo(const size_t ngridboxes,
                                    const Config &config);
/* return vector of dimensionless initial conditions
for thermodyanmic variables (p, temp, qv, qc) to
initialise cvode thermodynamics solver */

void run_cvodecoupld(const Config &config,
                     const RunSDMStep<auto, auto, auto> &sdm,
                     const CreateDetectorsPtr auto &dtrs,
                     const int t_end, const int couplstep)
/* create CVODE thermodynamics solver, superdroplets and gridboxes and
then run superdroplet model (SDM) coupled to the thermodynamics solver */
{
  /* CVODE thermodynamics solver */
  CvodeThermoSolver cvode(config,
                          initcvodethermo(sdm.gbxmaps.ngridboxes, config));

  /* vector containing all superdroplets within a
  struct that also holds their associated gridbox index.
  (all superdroplets have same solute properties) */
  const auto solute(std::make_shared<const SoluteProperties>());
  Kokkos::vector<SuperdropWithGbxindex> SDsInGBxs =
      create_superdrops_from_initSDsfile(config.initSDs_filename,
                                         config.nSDsvec,
                                         config.SDnspace, solute);

  /* vector containing all gridboxes in SDM domain */
  Kokkos::vector<GridBox> gridboxes =
      create_gridboxes(sdm.gbxmaps, dtrs, SDsInGBxs);

  /* prepare coupled model for timestepping */
  auto genpool = preparetotimestep(sdm,
                                   cvode, gridboxes,
                                   config.wetradiiinit,
                                   t_end, couplstep);

  /* run coupled model from t=0 to t=t_end */
  timestep_cvodecoupld(t_end, couplstep, sdm, cvode,
                    genpool, gridboxes, SDsInGBxs);

  std::cout << "\n ---- CVODE-SDM Coupled Run Complete ---- \n";
}

template <class MSDs, SdmProcess P, Observer O>
Kokkos::Random_XorShift64_Pool<>
preparetotimestep(const RunSDMStep<MSDs, P, O> &sdm,
                  CvodeThermoSolver &cvode,
                  Kokkos::vector<GridBox> &gridboxes,
                  const bool wetradiiinit,
                  const int t_end,
                  const int couplstep)
/* print some details about the cvode thermodynamics solver setup
and return pool of Kokkos' random number generator. Call
function to set superdroplet radii to equilibrium wet radius
if wetradiiinit is true. Call sdm observer's prepare function */
{
  const size_t ngbxs(gridboxes.size());

  cvode.print_init_ODEdata(step2dimlesstime(couplstep),
                           step2dimlesstime(t_end));

  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    set_thermostate(ii, cvode, gridboxes(ii).state);
  }

  if (wetradiiinit)
  {
    set_superdroplets_to_wetradius(gridboxes);
  }

  sdm.observer.prepare(sdm.logbooks);
  
  return Kokkos::Random_XorShift64_Pool<>(std::random_device{}()); // pool of Kokkos' random number generators
}

#endif // RUN_CVODECOUPLD_HPP