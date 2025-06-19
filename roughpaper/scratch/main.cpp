/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main.cpp
 * Project: scratch
 * Created Date: Monday 29th January 2024
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
 * rough paper for checking small things
 */

#include <Kokkos_Core.hpp>
#include <array>
#include <concepts>
#include <iostream>

#include "zarr/dataset.hpp"
#include "./cleotypes_sizes.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_movement.hpp"
#include "configuration/config.hpp"
#include "coupldyn_fromfile/fromfile_cartesian_dynamics.hpp"
#include "coupldyn_fromfile/fromfilecomms.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "initialise/init_all_supers_from_binary.hpp"
#include "initialise/initgbxsnull.hpp"
#include "initialise/initialconditions.hpp"
#include "initialise/timesteps.hpp"
#include "observers/gbxindex_observer.hpp"
#include "observers/massmoments_observer.hpp"
#include "observers/nsupers_observer.hpp"
#include "observers/observers.hpp"
#include "observers/state_observer.hpp"
#include "observers/streamout_observer.hpp"
#include "observers/superdrops_observer.hpp"
#include "observers/thermo_observer.hpp"
#include "observers/time_observer.hpp"
#include "observers/totnsupers_observer.hpp"
#include "observers/windvel_observer.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "zarr/fsstore.hpp"

template <typename Store>
inline Observer auto create_superdrops_observer(const Config &config, const Timesteps &tsteps,
                                                Dataset<Store> &dataset) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();

  CollectDataForDataset<Store> auto sdid = CollectSdId(dataset, maxchunk);
  CollectDataForDataset<Store> auto sdgbxindex = CollectSdgbxindex(dataset, maxchunk);
  CollectDataForDataset<Store> auto xi = CollectXi(dataset, maxchunk);
  CollectDataForDataset<Store> auto radius = CollectRadius(dataset, maxchunk);
  CollectDataForDataset<Store> auto msol = CollectMsol(dataset, maxchunk);
  CollectDataForDataset<Store> auto coord3 = CollectCoord3(dataset, maxchunk);
  CollectDataForDataset<Store> auto coord1 = CollectCoord1(dataset, maxchunk);
  CollectDataForDataset<Store> auto coord2 = CollectCoord2(dataset, maxchunk);

  const auto collect_data =
      coord1 >> coord2 >> coord3 >> msol >> radius >> xi >> sdgbxindex >> sdid;
  return SuperdropsObserver(obsstep, dataset, maxchunk, collect_data);
}

template <typename Store>
inline Observer auto create_gridbox_observer(const Config &config, const Timesteps &tsteps,
                                             Dataset<Store> &dataset) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();
  const auto ngbxs = config.get_ngbxs();

  const CollectDataForDataset<Store> auto thermo = CollectThermo(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Store> auto windvel = CollectWindVel(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Store> auto nsupers = CollectNsupers(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Store> auto collect_data = nsupers >> windvel >> thermo;
  return WriteToDatasetObserver(obsstep, dataset, collect_data);

  // const Observer auto obst = ThermoObserver(obsstep, dataset, maxchunk, ngbxs);
  // const Observer auto obsw = WindVelObserver(obsstep, dataset, maxchunk, ngbxs);
  // return obsw >> obst;

  // const Observer auto obsx = StateObserver(obsstep, dataset, maxchunk, ngbxs);
  // return obsx;

  // const Observer auto obsn = NsupersObserver(obsstep, dataset, maxchunk, ngbxs);
  // return obsn;
}

template <typename Store>
inline Observer auto create_observer2(const Config &config, const Timesteps &tsteps,
                                      Dataset<Store> &dataset) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();
  const auto ngbxs = config.get_ngbxs();

  const Observer auto obs0 = TimeObserver(obsstep, dataset, maxchunk, &step2dimlesstime);
  const Observer auto obs1 = GbxindexObserver(dataset, maxchunk, ngbxs);
  const Observer auto obs2 = MassMomentsObserver(obsstep, dataset, maxchunk, ngbxs);
  const Observer auto obs3 = MassMomentsRaindropsObserver(obsstep, dataset, maxchunk, ngbxs);
  const Observer auto obs4 = TotNsupersObserver(obsstep, dataset, maxchunk);
  const Observer auto obsx = create_gridbox_observer(config, tsteps, dataset);
  const Observer auto obssd = create_superdrops_observer(config, tsteps, dataset);

  return obssd >> obsx >> obs4 >> obs3 >> obs2 >> obs1 >> obs0;
}

/* ---------------------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------------------- */

template <typename Store>
inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     Dataset<Store> &dataset) {
  const auto obsstep = tsteps.get_obsstep();

  const Observer auto obs0 = StreamOutObserver(obsstep, &step2realtime);

  const Observer auto obs1 = create_observer2(config, tsteps, dataset);

  return obs0 >> obs1;
}

inline auto create_movement(const CartesianMaps &gbxmaps) {
  const Motion<CartesianMaps> auto motion = NullMotion{};
  const BoundaryConditions<CartesianMaps> auto boundary_conditions = NullBoundaryConditions{};

  return cartesian_movement(gbxmaps, motion, boundary_conditions);
}

template <GridboxMaps GbxMaps>
inline InitialConditions auto create_initconds(const Config &config, const GbxMaps &gbxmaps) {
  const auto initsupers = InitAllSupersFromBinary(config.get_initsupersfrombinary());
  const auto initgbxs = InitGbxsNull(gbxmaps.get_local_ngridboxes_hostcopy());

  return InitConds(initsupers, initgbxs);
}

inline CoupledDynamics auto create_coupldyn(const Config &config, const CartesianMaps &gbxmaps,
                                            const unsigned int couplstep,
                                            const unsigned int t_end) {
  const auto h_ndims = gbxmaps.get_global_ndims_hostcopy();
  const std::array<size_t, 3> ndims({h_ndims(0), h_ndims(1), h_ndims(2)});

  const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);

  return FromFileDynamics(config.get_fromfiledynamics(), couplstep, ndims, nsteps);
}

template <typename Store>
inline auto create_sdm(const Config &config, const Timesteps &tsteps, Dataset<Store> &dataset) {
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps = create_cartesian_maps(
      config.get_ngbxs(), config.get_nspacedims(), config.get_grid_filename());
  const MicrophysicalProcess auto microphys = NullMicrophysicalProcess{};
  const MoveSupersInDomain movesupers(create_movement(gbxmaps));
  const Observer auto obs = create_observer(config, tsteps, dataset);
  return SDMMethods(couplstep, gbxmaps, microphys, movesupers, obs);
}

int main(int argc, char *argv[]) {
  // print_type_sizes(argc, argv);

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
  const Config config("/home/m/m300950/CLEO/roughpaper/scratch/share/config.yaml");

  Kokkos::initialize(config.get_kokkos_initialization_settings());
  {
    Kokkos::print_configuration(std::cout);

    /* Create timestepping parameters from configuration */
    const Timesteps tsteps(config.get_timesteps());

    /* Create zarr store for writing output to storage */
    auto store = FSStore(config.get_zarrbasedir());
    auto dataset = Dataset(store);

    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm(create_sdm(config, tsteps, dataset));

    /* Solver of dynamics coupled to CLEO SDM */
    CoupledDynamics auto coupldyn(
        create_coupldyn(config, sdm.gbxmaps, tsteps.get_couplstep(), tsteps.get_t_end()));

    /* coupling between coupldyn and SDM */
    const CouplingComms<CartesianMaps, FromFileDynamics> auto comms = FromFileComms{};

    /* Initial conditions for CLEO run */
    const InitialConditions auto initconds = create_initconds(config, sdm.gbxmaps);

    /* Run CLEO (SDM coupled to dynamics solver) */
    const RunCLEO runcleo(sdm, coupldyn, comms);
    runcleo(initconds, tsteps.get_t_end());
  }
  Kokkos::finalize();

  const auto ttot = double{kokkostimer.seconds()};
  std::cout << "-------------------------------\n"
               "Total Program Duration: "
            << ttot << "s \n-------------------------------\n";

  MPI_Finalize();

  return 0;
}
/* ----------------------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------------*/
