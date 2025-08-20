/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_adia0d.cpp
 * Project: src
 * Created Date: Monday 29th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * runs the CLEO super-droplet model (SDM) for adiabatic parcel model example.
 * After make/compiling, execute for example via:
 * ./src/adia0d ../src/config/config.yaml
 */

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <stdexcept>
#include <string_view>

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_movement.hpp"
#include "configuration/communicator.hpp"
#include "configuration/config.hpp"
#include "coupldyn_cvode/cvodecomms.hpp"
#include "coupldyn_cvode/cvodedynamics.hpp"
#include "coupldyn_cvode/initgbxs_cvode.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "initialise/init_all_supers_from_binary.hpp"
#include "initialise/initialconditions.hpp"
#include "initialise/timesteps.hpp"
#include "observers/collect_data_for_simple_dataset.hpp"
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
#include "zarr/fsstore.hpp"
#include "zarr/simple_dataset.hpp"

inline CoupledDynamics auto create_coupldyn(const Config &config, const unsigned int couplstep) {
  return CvodeDynamics(config.get_cvodedynamics(), couplstep, &step2dimlesstime);
}

template <GridboxMaps GbxMaps>
inline InitialConditions auto create_initconds(const Config &config, const GbxMaps &gbxmaps) {
  const auto initsupers = InitAllSupersFromBinary(config.get_initsupersfrombinary());
  const auto initgbxs = InitGbxsCvode(config.get_cvodedynamics());

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
  const auto c = config.get_condensation();

  return Condensation(tsteps.get_condstep(), &step2dimlesstime, c.do_alter_thermo, c.maxniters,
                      c.rtol, c.atol, c.MINSUBTSTEP, &realtime2dimless);
}

template <typename Dataset, typename Store>
inline Observer auto create_superdrops_observer(const unsigned int interval, Dataset &dataset,
                                                Store &store, const int maxchunk) {
  CollectDataForDataset<Dataset> auto sdid = CollectSdId(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto sdgbxindex = CollectSdgbxindex(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto xi = CollectXi(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto radius = CollectRadius(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto msol = CollectMsol(dataset, maxchunk);

  const auto collect_sddata = msol >> radius >> xi >> sdgbxindex >> sdid;
  return SuperdropsObserver(interval, dataset, store, maxchunk, collect_sddata);
}

template <typename Dataset, typename Store>
inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     Dataset &dataset, Store &store) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();
  const auto ngbxs = config.get_ngbxs();

  const Observer auto obs1 = StreamOutObserver(obsstep * 10, &step2realtime);

  const Observer auto obs2 = TimeObserver(obsstep, dataset, store, maxchunk, &step2dimlesstime);

  const Observer auto obs3 = GbxindexObserver(dataset, store, maxchunk, ngbxs);

  const Observer auto obs4 = StateObserver(obsstep, dataset, maxchunk, ngbxs);

  const Observer auto obs5 = create_superdrops_observer(obsstep, dataset, store, maxchunk);

  return obs5 >> obs4 >> obs3 >> obs2 >> obs1;
}

template <typename Dataset, typename Store>
inline auto create_sdm(const Config &config, const Timesteps &tsteps, Dataset &dataset,
                       Store &store) {
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps = create_gbxmaps(config);
  const MicrophysicalProcess auto microphys = create_microphysics(config, tsteps);
  const MoveSupersInDomain movesupers = create_movement(gbxmaps);
  const Observer auto obs = create_observer(config, tsteps, dataset, store);

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

  /* Initialize Communicator here */
  init_communicator init_comm(argc, argv, config);

  /* Prevent this example from running with more than one MPI process */
  const auto comm_size = init_communicator::get_comm_size();
  if (comm_size > 1) {
    throw std::invalid_argument(
        "ERROR: The current example is not prepared to be run with more than one MPI process");
  }

  /* Initialise Kokkos device and host environments */
  Kokkos::initialize(config.get_kokkos_initialization_settings());
  {
    Kokkos::print_configuration(std::cout);

    /* Create timestepping parameters from configuration */
    const Timesteps tsteps(config.get_timesteps());

    /* Create Xarray dataset wit Zarr backend for writing output data to a store */
    auto store = FSStore(config.get_zarrbasedir());
    auto dataset = SimpleDataset(store);

    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm = create_sdm(config, tsteps, dataset, store);

    /* Create coupldyn solver and coupling between coupldyn and SDM */
    CoupledDynamics auto coupldyn = create_coupldyn(config, tsteps.get_couplstep());
    const CouplingComms<CartesianMaps, CvodeDynamics> auto comms = CvodeComms{};

    /* Initial conditions for CLEO run */
    const InitialConditions auto initconds = create_initconds(config, sdm.gbxmaps);

    /* Run CLEO (SDM coupled to dynamics solver) */
    const RunCLEO runcleo(sdm, coupldyn, comms);
    runcleo(initconds, tsteps.get_t_end());
  }
  Kokkos::finalize();

  const auto ttot = double{kokkostimer.seconds()};
  std::cout << "-----\n Total Program Duration: " << ttot << "s \n-----\n";

  return 0;
}
