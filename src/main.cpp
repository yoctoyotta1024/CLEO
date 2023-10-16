/*
 * ----- CLEO -----
 * File: main.cpp
 * Project: src
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 16th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * runs the CLEO super-droplet model (SDM)
 * after make/compiling, execute for example via:
 * ./src/runCLEO ../src/config/config.txt
 */

#include <iostream>
#include <stdexcept>
#include <string_view>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "coupldyn_fromfile/fromfiledynamics.hpp"
#include "initialise/config.hpp"
#include "initialise/timesteps.hpp"
#include "observers/constintervalobs.hpp"
#include "observers/observers.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"
#include "sdmdomain/gridboxmaps.hpp"
#include "zarr/fsstore.hpp"

CoupledDynamics auto
create_coupldyn(const Config &config,
                const unsigned int coupldynstep)
{
  return FromFileDynamics(config, coupldynstep); 
}

Observer auto
create_observer(const unsigned int obsstep)
{
  return ConstIntervalObs(obsstep); 
}

auto create_sdm(const Config &config,
                const Timesteps &tsteps,
                const CoupledDynamics auto &coupldyn)
{
  const GridboxMaps gbxmaps(config);
  const MicrophysicsProcess microphys;
  const MoveSupersInDomain movesupers(tsteps.get_motionstep());
  const Observer auto obs(create_observer(tsteps.get_obsstep()));

  return SDMMethods(coupldyn, gbxmaps,
                    microphys, movesupers, obs);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    throw std::invalid_argument("configuration file(s) not specified");
  }
 
  Kokkos::Timer kokkostimer;

  /* Read input parameters from configuration file(s) */
  const std::string_view config_filename(argv[1]);    // path to configuration file
  const Config config(config_filename);
  const Timesteps tsteps(config); // timesteps for model (e.g. coupling and end time)

  /* Create zarr store for writing output to storage */
  FSStore fsstore(config.zarrbasedir);
    
  /* Solver of dynamics coupled to CLEO SDM */
  const CoupledDynamics auto coupldyn(
      create_coupldyn(config, tsteps.get_couplstep()));

  /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
  const SDMMethods sdm(create_sdm(config, tsteps, coupldyn));

  /* Run CLEO (SDM coupled to dynamics solver) */
  Kokkos::initialize(argc, argv);
  {
    const RunCLEO runcleo(sdm, coupldyn);
    runcleo(tsteps.get_t_end());
  }
  Kokkos::finalize();
  
  // const double ttot(kokkostimer.seconds());
  // std::cout << "-----\n Total Program Duration: "
  //           << ttot << "s \n-----\n";

  return 0;
}
