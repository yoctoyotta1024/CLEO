/*
 * ----- CLEO -----
 * File: main.cpp
 * Project: src
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
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

#include <Kokkos_Core.hpp>

#include "initialise/config.hpp"
#include "initialise/timesteps.hpp"
#include "zarr/fsstore.hpp"
#include "runcleo/runcleo.hpp"

CoupledDynamics
create_coupldyn(const Config &config,
                const unsigned int coupldynstep)
{
  return CoupledDynamics(config, coupldynstep); 
}

CLEOSDM
create_sdm(const Config &config,
           const Timesteps &tsteps,
           const unsigned int couplstep)
{
  GridboxMaps gbxmaps(config);
  MicrophysicsProcess microphys;
  MoveSupersInDomain movesupers(tsteps.get_motionstep());
  Observer obs(tsteps.get_obsstep());

  return CLEOSDM(gbxmaps, microphys,
                 movesupers, obs, couplstep);
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
  const auto coupldyn(create_coupldyn(config, tsteps.get_couplstep()));

  /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
  const auto sdm(create_sdm(config, tsteps, coupldyn.get_couplstep()));

  /* Run CLEO (SDM coupled to dynamics solver) */
  Kokkos::initialize(argc, argv);
  {
    run_cleo(tsteps.get_t_end(), sdm, coupldyn);
  }

  Kokkos::finalize();
  const double ttot(kokkostimer.seconds());
  std::cout << "-----\n Total Program Duration: "
            << ttot << "s \n-----\n";

  return 0;
}
