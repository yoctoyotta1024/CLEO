/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: main_adia0D.cpp
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
 * File Description:
 * runs the CLEO super-droplet model (SDM)
 * after make/compiling, execute for example via:
 * ./src/adia0D ../src/config/config.txt
 */

#include <concepts>
#include <iostream>
#include <string_view>
#include <stdexcept>

#include <Kokkos_Core.hpp>

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/cartesianmotion.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "coupldyn_cvode/cvodecomms.hpp"
#include "coupldyn_cvode/cvodedynamics.hpp"
#include "coupldyn_cvode/initgbxs_cvode.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "initialise/config.hpp"
#include "initialise/initsupers_frombinary.hpp"
#include "initialise/timesteps.hpp"
#include "observers/observers.hpp"
#include "observers/printobs.hpp"
#include "observers/stateobs.hpp"
#include "observers/supersattrsobs.hpp"
#include "observers/timeobs.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/initialconditions.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/condensation.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "zarr/fsstore.hpp"
#include "zarr/superdropattrsbuffers.hpp"
#include "zarr/superdropsbuffers.hpp"

inline CoupledDynamics auto create_coupldyn(const Config &config, const unsigned int couplstep) {
  return CvodeDynamics(config, couplstep, &step2dimlesstime);
}

inline InitialConditions auto create_initconds(const Config &config) {
  const InitSupersFromBinary initsupers(config);
  const InitGbxsCvode initgbxs(config);

  return InitConds(initsupers, initgbxs);
}

inline GridboxMaps auto create_gbxmaps(const Config &config) {
  const auto gbxmaps = create_cartesian_maps(config.ngbxs, config.nspacedims, config.grid_filename);
  return gbxmaps;
}

inline MicrophysicalProcess auto create_microphysics(const Config &config,
                                                     const Timesteps &tsteps) {
  const MicrophysicalProcess auto cond = Condensation(
      tsteps.get_condstep(), config.doAlterThermo, config.cond_iters, &step2dimlesstime,
      config.cond_rtol, config.cond_atol, config.cond_SUBTSTEP, &realtime2dimless);
  return cond;
}

inline Observer auto create_supersattrs_observer(const unsigned int interval, FSStore &store,
                                                 const int maxchunk) {
  SuperdropsBuffers auto buffers =
      SdIdBuffer() >> XiBuffer() >> MsolBuffer() >> RadiusBuffer() >> SdgbxindexBuffer();
  return SupersAttrsObserver(interval, store, maxchunk, buffers);
}

inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     FSStore &store) {
  const auto obsstep = (unsigned int)tsteps.get_obsstep();
  const auto maxchunk = int{config.maxchunk};

  const Observer auto obs1 = PrintObserver(obsstep * 100, &step2realtime);

  const Observer auto obs2 = TimeObserver(obsstep, store, maxchunk, &step2dimlesstime);

  const Observer auto obs3 = StateObserver(obsstep, store, maxchunk, config.ngbxs);

  const Observer auto obs4 = create_supersattrs_observer(obsstep, store, maxchunk);

  return obs1 >> obs2 >> obs3 >> obs4;
}

inline auto create_sdm(const Config &config, const Timesteps &tsteps, FSStore &store) {
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps(create_gbxmaps(config));
  const MicrophysicalProcess auto microphys(create_microphysics(config, tsteps));
  const Motion<CartesianMaps> auto movesupers = NullMotion{};
  const Observer auto obs(create_observer(config, tsteps, store));

  return SDMMethods(couplstep, gbxmaps, microphys, movesupers, obs);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    throw std::invalid_argument("configuration file(s) not specified");
  }

  Kokkos::Timer kokkostimer;

  /* Read input parameters from configuration file(s) */
  const std::string_view config_filename(argv[1]);  // path to configuration file
  const Config config(config_filename);
  const Timesteps tsteps(config);  // timesteps for model (e.g. coupling and end time)

  /* Create zarr store for writing output to storage */
  FSStore fsstore(config.zarrbasedir);

  /* create coupldyn solver and coupling between coupldyn and SDM */
  CoupledDynamics auto coupldyn(create_coupldyn(config, tsteps.get_couplstep()));
  const CouplingComms<CvodeDynamics> auto comms = CvodeComms{};

  /* Initial conditions for CLEO run */
  const InitialConditions auto initconds = create_initconds(config);

  /* Initialise Kokkos device and host environments */
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
  std::cout << "-----\n Total Program Duration: " << ttot << "s \n-----\n";

  return 0;
}
