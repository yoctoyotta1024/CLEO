/*
 * ----- CLEO -----
 * File: main_longcolls.cpp
 * Project: src
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 11th January 2024
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
 * ./src/longcolls ../src/config/config.txt
 */

#include <iostream>
#include <stdexcept>
#include <string_view>
#include <concepts>
#include <cmath>

#include <Kokkos_Core.hpp>

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/cartesianmotion.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"

#include "coupldyn_null/nulldynamics.hpp"
#include "coupldyn_null/nulldyncomms.hpp"

#include "gridboxes/gridboxmaps.hpp"

#include "initialise/config.hpp"
#include "initialise/timesteps.hpp"
#include "initialise/initgbxs_null.hpp"
#include "initialise/initsupers_frombinary.hpp"

#include "observers/observers.hpp"
#include "observers/printobs.hpp"
#include "observers/timeobs.hpp"
#include "observers/supersattrsobs.hpp"

#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/initialconditions.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"

#include "superdrops/coalescence.hpp"
#include "superdrops/collisionprobs/longhydroprob.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/microphysicalprocess.hpp"

#include "zarr/fsstore.hpp"
#include "zarr/superdropattrsbuffers.hpp"
#include "zarr/superdropsbuffers.hpp"

inline InitialConditions auto
create_initconds(const Config &config)
{
  const InitSupersFromBinary initsupers(config);
  const InitGbxsNull initgbxs(config);

  return InitConds(initsupers, initgbxs);
}

inline GridboxMaps auto
create_gbxmaps(const Config &config)
{
  const auto gbxmaps = create_cartesian_maps(config.ngbxs,
                                             config.nspacedims,
                                             config.grid_filename);
  return gbxmaps;
}

inline MicrophysicalProcess auto
create_microphysics(const Config &config, const Timesteps &tsteps)
{
  const PairProbability auto prob = LongHydroProb();
  const MicrophysicalProcess auto colls = CollCoal(tsteps.get_collstep(),
                                                   &step2realtime,
                                                   prob);                                                    
  return colls;
}

inline Motion<CartesianMaps> auto
create_motion(const unsigned int motionstep)
{
  return NullMotion{};                                                                               
}

inline Observer auto
create_supersattrs_observer(const unsigned int interval,
                            FSStore &store,
                            const int maxchunk)
{
  SuperdropsBuffers auto buffers = SdIdBuffer() >>
                                   XiBuffer() >>
                                   MsolBuffer() >>
                                   RadiusBuffer();
  return SupersAttrsObserver(interval, store, maxchunk, buffers);
}

inline Observer auto
create_observer(const Config &config,
                const Timesteps &tsteps,
                FSStore &store)
{
  const auto obsstep = (unsigned int)tsteps.get_obsstep();
  const auto maxchunk = int{config.maxchunk};

  const Observer auto obs1 = PrintObserver(obsstep, &step2realtime);

  const Observer auto obs2 = TimeObserver(obsstep, store, maxchunk,
                                          &step2dimlesstime);

  const Observer auto obs3 = create_supersattrs_observer(obsstep, store,
                                                          maxchunk);

  return obs1 >> obs2 >> obs3;
}

inline auto create_sdm(const Config &config,
                       const Timesteps &tsteps,
                       FSStore &store)
{
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps(create_gbxmaps(config));
  const MicrophysicalProcess auto microphys(create_microphysics(config, tsteps));
  const Motion<CartesianMaps> auto movesupers(create_motion(tsteps.get_motionstep()));
  const Observer auto obs(create_observer(config, tsteps, store));

  return SDMMethods(couplstep, gbxmaps,
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
  const std::string_view config_filename(argv[1]); // path to configuration file
  const Config config(config_filename);
  const Timesteps tsteps(config); // timesteps for model (e.g. coupling and end time)

  /* Create zarr store for writing output to storage */
  FSStore fsstore(config.zarrbasedir);

  /* create coupldyn solver and coupling between coupldyn and SDM */
  const CoupledDynamics auto
      coupldyn = NullDynamics(tsteps.get_couplstep());
  const CouplingComms<NullDynamics> auto comms = NullDynComms{};

  /* Initial conditions for CLEO run */
  const InitialConditions auto initconds = create_initconds(config);

  /* Initialise Kokkos parallel environment */
  Kokkos::initialize(argc, argv);
  {
    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm(create_sdm(config, tsteps, fsstore));

    /* Run CLEO (SDM coupled to dynamics solver) */
    const RunCLEO runcleo(sdm, coupldyn, comms);
    runcleo(initconds, tsteps.get_t_end());
  }
  Kokkos::finalize();

  const auto ttot = double{kokkostimer.seconds()};
  std::cout << "-----\n Total Program Duration: "
            << ttot << "s \n-----\n";

  return 0;
}