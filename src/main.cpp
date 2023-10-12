/*
 * ----- CLEO -----
 * File: main.cpp
 * Project: src
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 12th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * runs the CLEO super-droplet model (SDM)
 * after make/compiling, execute for example via: ./src/runCLEO "../src/config/config.txt"
 */

#include <iostream>
#include <stdexcept>
#include <string_view>

#include <Kokkos_Core.hpp>

struct Config
{
  Config(const std::string_view configfile) {}
};

struct FSStore
{
};

struct ThermoSolver
{
  int run() const
  {
    std::cout << "Hi World\n";

    return 0;
  }
};

struct GridboxMaps 
{
};

struct Microphys
{
};

struct Motion
{
};

struct Observer
{
};

class CLEOSDM
{
private:
  ThermoSolver solver;
  GridboxMaps gbxmaps; // maps from gridbox indexes to domain coordinates
  Microphys microphys; // microphysical process
  Motion motion; // super-droplets' motion
  Observer obs; // observer

public:
  int run() const
  {
    std::cout << "Hello World\n";

    return 0;
  }
};

int run_cleo(const CLEOSDM &sdm, const ThermoSolver &solver)
{
  solver.run();
  sdm.run();

  return 0;
}

int main(int argc, char *argv[])
{
  Kokkos::Timer kokkostimer;

  if (argc < 2)
  {
    throw std::invalid_argument("configuration file(s) not specified");
  }

  /* read input parameters from configuration file(s) */
  const std::string_view configfile = argv[1];    // path to configuration file
  const Config config(configfile);

  /* create zarr store for writing output to storage */
  FSStore fsstore;
  
  /* Thermodynamics Solver to couple to CLEO */
  const ThermoSolver solver;

  /* CLEO SDM coupled to a Thermodyanmics Solver */
  const CLEOSDM sdm;

  /* run coupled CLEO SDM with Kokkos */
  Kokkos::initialize(argc, argv);
  {
    run_cleo(sdm, solver);
  }

  Kokkos::finalize();
  std::cout << "  ------ Total Duration: " << kokkostimer.seconds() << "s ----- \n";

  return 0;
}
