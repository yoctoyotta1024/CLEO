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

struct TimeSteps
{
private:
  unsigned int t_end = 5;
public:
  unsigned int get_t_end() const { return t_end; }
};

struct GridBoxMaps 
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

struct GridBoxes
{
};

struct SuperDrops
{
};

struct CLEOSDM
{
  TimeSteps tsteps; // timesteps for model (e.g. coupling and end time)
  GridBoxMaps gbxmaps; // maps from gridbox indexes to domain coordinates
  Microphys microphys; // microphysical process
  Motion motion; // super-droplets' motion
  Observer obs; // observer

  GridBoxes generate_gridboxes() const
  {
    return GridBoxes{};
  }

  SuperDrops generate_superdrops() const
  {
    return SuperDrops{};
  }

  int prepare_to_timestep(const GridBoxes &GBxs,
                          const SuperDrops &SDs) const
  {
    return 0;
  }

  void run_step(const unsigned int t_mdl) const
  {
    std::cout << "SDM Call @ t=" << t_mdl << "\n" ;
  }
};

struct CoupledDynamics
{
  int prepare_to_timestep(const GridBoxes &GBxs,
                          const SuperDrops &SDs) const
  {
    return 0;
  }

  void run_step(const unsigned int t_mdl) const
  {
    std::cout << "Dyn Call @ t=" << t_mdl << "\n" ;
  }
};

struct RunStats
{
private:
  Kokkos::Timer kokkostimer;
  double time1;
  double time2;

public:
  void pre_timestepping()
  {
    time1 = kokkostimer.seconds();
  }

  void post_timestepping()
  {
    time2 = kokkostimer.seconds();
  }

  void summary()
  {
    std::cout << "\n ----- CLEO run complete ----- \n"
              << "       Duration: " << time2 << "s \n"
              << "       Initialisation: " << time1 << "s \n"
              << "       Timestepping: " << time2 - time1 << "s \n"
              << "------------------------------- \n";
  }
};

void start_step(const unsigned int t_mdl)
{}

unsigned int next_step(unsigned int t_mdl)
{
  return ++t_mdl;
}

int run_cleo(const CLEOSDM &sdm, const CoupledDynamics &coupldyn)
{
  // generate runtime objects
  RunStats stats;
  GridBoxes GBxs(sdm.generate_gridboxes());
  SuperDrops SDs(sdm.generate_superdrops());

  // prepare CLEO for timestepping
  coupldyn.prepare_to_timestep(GBxs, SDs);
  sdm.prepare_to_timestep(GBxs, SDs);
  stats.pre_timestepping();
  
  // timestep CLEO from t=0 to t=t_end
  int t_mdl = 0;
  while (t_mdl <= sdm.tsteps.get_t_end())
  {
    /* start step (in general involves coupling) */
    start_step(t_mdl);

    /* advance SDM (optionally concurrent to dynamics solver) */
    sdm.run_step(t_mdl);

    /* advance dynamics solver (optionally concurrent to SDM) */
    coupldyn.run_step(t_mdl);

    /* proceed to next step (in general involves coupling) */
    t_mdl = next_step(t_mdl);
  }
  stats.post_timestepping();
  
  // summary of runtime statistics
  stats.summary();

  return 0;
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    throw std::invalid_argument("configuration file(s) not specified");
  }
 
  Kokkos::Timer kokkostimer;

  /* read input parameters from configuration file(s) */
  const std::string_view configfile = argv[1];    // path to configuration file
  const Config config(configfile);

  /* create zarr store for writing output to storage */
  FSStore fsstore;
  
  /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
  const CLEOSDM sdm;

  /* Solver of Dynamics coupled to CLEO SDM */
  const CoupledDynamics coupldyn;

  /* run coupled CLEO SDM with Kokkos */
  Kokkos::initialize(argc, argv);
  {
    run_cleo(sdm, coupldyn);
  }

  Kokkos::finalize();
  const double ttot(kokkostimer.seconds());
  std::cout << "-----\n Total Program Duration: "
            << ttot << "s \n-----\n";

  return 0;
}
