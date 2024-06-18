/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_adia0D.cpp
 * Project: src
 * Created Date: Monday 29th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 18th June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * runs the CLEO super-droplet model (SDM) for adiabatic parcel model example.
 * After make/compiling, execute for example via:
 * ./src/adia0D ../src/config/config.yaml
 */

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <stdexcept>
#include <string_view>

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/cartesianmotion.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "cartesiandomain/null_boundary_conditions.hpp"
#include "coupldyn_cvode/cvodecomms.hpp"
#include "coupldyn_cvode/cvodedynamics.hpp"
#include "coupldyn_cvode/initgbxs_cvode.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "initialise/config.hpp"
#include "initialise/init_all_supers_from_binary.hpp"
#include "initialise/initialconditions.hpp"
#include "initialise/timesteps.hpp"
#include "observers/gbxindex_observer.hpp"
#include "observers/observers.hpp"
#include "observers/state_observer.hpp"
#include "observers/streamout_observer.hpp"
#include "observers/superdrops_observer.hpp"
#include "observers/time_observer.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/condensation.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "zarr/dataset.hpp"
#include "zarr/fsstore.hpp"

inline CoupledDynamics auto create_coupldyn(const Config &config, const unsigned int couplstep) {
  return CvodeDynamics(config.get_cvodedynamics(), couplstep, &step2dimlesstime);
}

inline InitialConditions auto create_initconds(const Config &config) {
  const InitAllSupersFromBinary initsupers(config.get_initsupersfrombinary());
  const InitGbxsCvode initgbxs(config.get_cvodedynamics());

  return InitConds(initsupers, initgbxs);
}

inline GridboxMaps auto create_gbxmaps(const Config &config) {
  const auto gbxmaps = create_cartesian_maps(config.get_ngbxs(), config.get_nspacedims(),
                                             config.get_grid_filename());
  return gbxmaps;
}

inline auto create_movement(const CartesianMaps &gbxmaps) {
  const Motion<CartesianMaps> auto motion = NullMotion{};
  const auto boundary_conditions = NullBoundaryConditions{};

  return MoveSupersInDomain(gbxmaps, motion, boundary_conditions);
}

inline MicrophysicalProcess auto create_microphysics(const Config &config,
                                                     const Timesteps &tsteps) {
  const auto c = config.get_condensation();

  return Condensation(tsteps.get_condstep(), &step2dimlesstime, c.do_alter_thermo, c.maxniters,
                      c.rtol, c.atol, c.MINSUBTSTEP, &realtime2dimless);
}

template <typename Store>
inline Observer auto create_superdrops_observer(const unsigned int interval,
                                                Dataset<Store> &dataset, const int maxchunk) {
  CollectDataForDataset<Store> auto sdid = CollectSdId(dataset, maxchunk);
  CollectDataForDataset<Store> auto sdgbxindex = CollectSdgbxindex(dataset, maxchunk);
  CollectDataForDataset<Store> auto xi = CollectXi(dataset, maxchunk);
  CollectDataForDataset<Store> auto radius = CollectRadius(dataset, maxchunk);
  CollectDataForDataset<Store> auto msol = CollectMsol(dataset, maxchunk);

  const auto collect_sddata = msol >> radius >> xi >> sdgbxindex >> sdid;
  return SuperdropsObserver(interval, dataset, maxchunk, collect_sddata);
}

template <typename Store>
inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     Dataset<Store> &dataset) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();
  const auto ngbxs = config.get_ngbxs();

  const Observer auto obs1 = StreamOutObserver(obsstep * 10, &step2realtime);

  const Observer auto obs2 = TimeObserver(obsstep, dataset, maxchunk, &step2dimlesstime);

  const Observer auto obs3 = GbxindexObserver(dataset, maxchunk, ngbxs);

  const Observer auto obs4 = StateObserver(obsstep, dataset, maxchunk, ngbxs);

  const Observer auto obs5 = create_superdrops_observer(obsstep, dataset, maxchunk);

  return obs5 >> obs4 >> obs3 >> obs2 >> obs1;
}

template <typename Store>
inline auto create_sdm(const Config &config, const Timesteps &tsteps, Dataset<Store> &dataset) {
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps(create_gbxmaps(config));
  const MicrophysicalProcess auto microphys(create_microphysics(config, tsteps));
  const MoveSupersInDomain movesupers(create_movement(gbxmaps));
  const Observer auto obs(create_observer(config, tsteps, dataset));

  return SDMMethods(couplstep, gbxmaps, microphys, movesupers, obs);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    throw std::invalid_argument("configuration file(s) not specified");
  }

  Kokkos::Timer kokkostimer;

  /* Read input parameters from configuration file(s) */
  const std::filesystem::path config_filename(argv[1]);  // path to configuration file
  const Config config(config_filename);
  const Timesteps tsteps(config.get_timesteps());

  /* Create Xarray dataset wit Zarr backend for writing output data to a store */
  auto store = FSStore(config.get_zarrbasedir());
  auto dataset = Dataset(store);

  /* Create coupldyn solver and coupling between coupldyn and SDM */
  CoupledDynamics auto coupldyn(create_coupldyn(config, tsteps.get_couplstep()));
  const CouplingComms<CvodeDynamics> auto comms = CvodeComms{};

  /* Initial conditions for CLEO run */
  const InitialConditions auto initconds = create_initconds(config);

  /* Initialise Kokkos device and host environments */
  Kokkos::initialize(argc, argv);
  {
    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm(create_sdm(config, tsteps, dataset));

    /* Run CLEO (SDM coupled to dynamics solver) */
    const RunCLEO runcleo(sdm, coupldyn, comms);
    runcleo(initconds, tsteps.get_t_end());
  }
  Kokkos::finalize();

  const auto ttot = double{kokkostimer.seconds()};
  std::cout << "-----\n Total Program Duration: " << ttot << "s \n-----\n";

  return 0;
}
