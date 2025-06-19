/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_golcolls.cpp
 * Project: src
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 28th May 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * runs the CLEO super-droplet model (SDM) for 0-D box model with Golovin's kernel.
 * After make/compiling, execute for example via:
 * ./src/golcolls ../src/config/config.yaml
 */

#include <Kokkos_Core.hpp>
#include <cmath>
#include <concepts>
#include <iostream>
#include <stdexcept>
#include <string_view>

#include "zarr/dataset.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_movement.hpp"
#include "configuration/config.hpp"
#include "coupldyn_null/nulldynamics.hpp"
#include "coupldyn_null/nulldyncomms.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "initialise/init_all_supers_from_binary.hpp"
#include "initialise/initgbxsnull.hpp"
#include "initialise/initialconditions.hpp"
#include "initialise/timesteps.hpp"
#include "observers/gbxindex_observer.hpp"
#include "observers/observers.hpp"
#include "observers/streamout_observer.hpp"
#include "observers/superdrops_observer.hpp"
#include "observers/time_observer.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/golovinprob.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "zarr/fsstore.hpp"

template <GridboxMaps GbxMaps>
inline InitialConditions auto create_initconds(const Config &config, const GbxMaps &gbxmaps) {
  const auto initsupers = InitAllSupersFromBinary(config.get_initsupersfrombinary());
  const auto initgbxs = InitGbxsNull(gbxmaps.get_local_ngridboxes_hostcopy());

  return InitConds(initsupers, initgbxs);
}

inline GridboxMaps auto create_gbxmaps(const Config &config) {
  const auto gbxmaps = create_cartesian_maps(config.get_ngbxs(), config.get_nspacedims(),
                                             config.get_grid_filename());
  return gbxmaps;
}

inline auto create_movement(const CartesianMaps &gbxmaps) {
  const Motion<CartesianMaps> auto motion = NullMotion{};
  const BoundaryConditions<CartesianMaps> auto boundary_conditions = NullBoundaryConditions{};

  return cartesian_movement(gbxmaps, motion, boundary_conditions);
}

inline MicrophysicalProcess auto create_microphysics(const Config &config,
                                                     const Timesteps &tsteps) {
  const PairProbability auto prob = GolovinProb();
  const MicrophysicalProcess auto colls = CollCoal(tsteps.get_collstep(), &step2realtime, prob);
  return colls;
}

template <typename Store>
inline Observer auto create_superdrops_observer(const unsigned int interval,
                                                Dataset<Store> &dataset, const int maxchunk) {
  CollectDataForDataset<Store> auto sdid = CollectSdId(dataset, maxchunk);
  CollectDataForDataset<Store> auto xi = CollectXi(dataset, maxchunk);
  CollectDataForDataset<Store> auto radius = CollectRadius(dataset, maxchunk);
  CollectDataForDataset<Store> auto msol = CollectMsol(dataset, maxchunk);

  const auto collect_sddata = msol >> radius >> xi >> sdid;
  return SuperdropsObserver(interval, dataset, maxchunk, collect_sddata);
}

template <typename Store>
inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     Dataset<Store> &dataset) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();

  const Observer auto obs0 = StreamOutObserver(obsstep, &step2realtime);

  const Observer auto obs1 = TimeObserver(obsstep, dataset, maxchunk, &step2dimlesstime);

  const Observer auto obssd = create_superdrops_observer(obsstep, dataset, maxchunk);

  return obssd >> obs1 >> obs0;
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

  MPI_Init(&argc, &argv);

  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size > 1) {
    std::cout << "ERROR: The current example is not prepared"
              << " to be run with more than one MPI process" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  Kokkos::Timer kokkostimer;

  /* Read input parameters from configuration file(s) */
  const std::filesystem::path config_filename(argv[1]);  // path to configuration file
  const Config config(config_filename);

  /* Initialise Kokkos parallel environment */
  Kokkos::initialize(config.get_kokkos_initialization_settings());
  {
    Kokkos::print_configuration(std::cout);

    /* Create timestepping parameters from configuration */
    const Timesteps tsteps(config.get_timesteps());

    /* Create Xarray dataset wit Zarr backend for writing output data to a store */
    auto store = FSStore(config.get_zarrbasedir());
    auto dataset = Dataset(store);

    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm(create_sdm(config, tsteps, dataset));

    /* Create coupldyn solver and coupling between coupldyn and SDM */
    const CoupledDynamics auto coupldyn = NullDynamics(tsteps.get_couplstep());
    const CouplingComms<CartesianMaps, NullDynamics> auto comms = NullDynComms{};

    /* Initial conditions for CLEO run */
    const InitialConditions auto initconds = create_initconds(config, sdm.gbxmaps);

    /* Run CLEO (SDM coupled to dynamics solver) */
    const RunCLEO runcleo(sdm, coupldyn, comms);
    runcleo(initconds, tsteps.get_t_end());
  }
  Kokkos::finalize();

  const auto ttot = double{kokkostimer.seconds()};
  std::cout << "-----\n Total Program Duration: " << ttot << "s \n-----\n";

  MPI_Finalize();

  return 0;
}
