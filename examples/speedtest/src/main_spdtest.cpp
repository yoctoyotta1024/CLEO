/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_spdtest.cpp
 * Project: src
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * runs the CLEO super-droplet model (SDM) for speed test model example.
 * After make/compiling, execute for example via:
 * ./src/spdtest ../src/config/config.yaml
 */

#include <Kokkos_Core.hpp>
#include <cmath>
#include <concepts>
#include <iostream>
#include <stdexcept>
#include <string_view>

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_motion.hpp"
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
#include "observers/collect_data_for_simple_dataset.hpp"
#include "observers/gbxindex_observer.hpp"
#include "observers/massmoments_observer.hpp"
#include "observers/nsupers_observer.hpp"
#include "observers/observers.hpp"
#include "observers/state_observer.hpp"
#include "observers/streamout_observer.hpp"
#include "observers/superdrops_observer.hpp"
#include "observers/time_observer.hpp"
#include "observers/totnsupers_observer.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/longhydroprob.hpp"
#include "superdrops/condensation.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/terminalvelocity.hpp"
#include "zarr/fsstore.hpp"
#include "zarr/simple_dataset.hpp"

inline CoupledDynamics auto create_coupldyn(const Config &config, const CartesianMaps &gbxmaps,
                                            const unsigned int couplstep,
                                            const unsigned int t_end) {
  const auto h_ndims = gbxmaps.get_global_ndims_hostcopy();
  const std::array<size_t, 3> ndims({h_ndims(0), h_ndims(1), h_ndims(2)});

  const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);

  return FromFileDynamics(config.get_fromfiledynamics(), couplstep, ndims, nsteps);
}

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

inline auto create_movement(const unsigned int motionstep, const CartesianMaps &gbxmaps) {
  const auto terminalv = RogersGKTerminalVelocity{};
  const Motion<CartesianMaps> auto motion =
      CartesianMotion(motionstep, &step2dimlesstime, terminalv);

  const BoundaryConditions<CartesianMaps> auto boundary_conditions = NullBoundaryConditions{};

  return cartesian_movement(gbxmaps, motion, boundary_conditions);
}

inline MicrophysicalProcess auto config_condensation(const Config &config,
                                                     const Timesteps &tsteps) {
  const auto c = config.get_condensation();

  return Condensation(tsteps.get_condstep(), &step2dimlesstime, c.do_alter_thermo, c.maxniters,
                      c.rtol, c.atol, c.MINSUBTSTEP, &realtime2dimless);
}

inline MicrophysicalProcess auto config_collisions(const Config &config, const Timesteps &tsteps) {
  const PairProbability auto coalprob = LongHydroProb(1.0);
  const MicrophysicalProcess auto coal = CollCoal(tsteps.get_collstep(), &step2realtime, coalprob);

  return coal;
}

inline MicrophysicalProcess auto create_microphysics(const Config &config,
                                                     const Timesteps &tsteps) {
  const MicrophysicalProcess auto cond = config_condensation(config, tsteps);

  const MicrophysicalProcess auto colls = config_collisions(config, tsteps);

  return colls >> cond;
}

template <typename Dataset, typename Store>
inline Observer auto create_superdrops_observer(const unsigned int interval, Dataset &dataset,
                                                Store &store, const int maxchunk) {
  CollectDataForDataset<Dataset> auto sdid = CollectSdId(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto sdgbxindex = CollectSdgbxindex(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto xi = CollectXi(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto radius = CollectRadius(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto msol = CollectMsol(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto coord3 = CollectCoord3(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto coord1 = CollectCoord1(dataset, maxchunk);
  CollectDataForDataset<Dataset> auto coord2 = CollectCoord2(dataset, maxchunk);

  const auto collect_sddata =
      coord2 >> coord1 >> coord3 >> msol >> radius >> xi >> sdgbxindex >> sdid;
  return SuperdropsObserver(interval, dataset, store, maxchunk, collect_sddata);
}

template <typename Dataset>
inline Observer auto create_gridboxes_observer(const unsigned int interval, Dataset &dataset,
                                               const int maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Dataset> auto thermo = CollectThermo(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Dataset> auto windvel = CollectWindVel(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Dataset> auto nsupers = CollectNsupers(dataset, maxchunk, ngbxs);

  const CollectDataForDataset<Dataset> auto collect_gbxdata = nsupers >> windvel >> thermo;
  return WriteToDatasetObserver(interval, dataset, collect_gbxdata);
}

template <typename Dataset, typename Store>
inline Observer auto create_bulk_observer(const unsigned int interval, Dataset &dataset,
                                          Store &store, const int maxchunk, const size_t ngbxs) {
  const Observer auto obs2 = TimeObserver(interval, dataset, store, maxchunk, &step2dimlesstime);

  const Observer auto obs3 = GbxindexObserver(dataset, store, maxchunk, ngbxs);

  const Observer auto obs4 = TotNsupersObserver(interval, dataset, store, maxchunk);

  const Observer auto obs5 = MassMomentsObserver(interval, dataset, store, maxchunk, ngbxs);

  const Observer auto obs6 =
      MassMomentsRaindropsObserver(interval, dataset, store, maxchunk, ngbxs);

  const Observer auto obsgbx = create_gridboxes_observer(interval, dataset, maxchunk, ngbxs);

  return obsgbx >> obs6 >> obs5 >> obs4 >> obs3 >> obs2;
}

template <typename Dataset, typename Store>
inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     Dataset &dataset, Store &store) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();

  const Observer auto obs0 = StreamOutObserver(obsstep * 10, &step2realtime);

  const Observer auto obsblk =
      create_bulk_observer(obsstep, dataset, store, maxchunk, config.get_ngbxs());

  const Observer auto obssd = create_superdrops_observer(obsstep, dataset, store, maxchunk);

  return obssd >> obsblk >> obs0;
}

template <typename Dataset, typename Store>
inline auto create_sdm(const Config &config, const Timesteps &tsteps, Dataset &dataset,
                       Store &store) {
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps(create_gbxmaps(config));
  const MicrophysicalProcess auto microphys(create_microphysics(config, tsteps));
  const MoveSupersInDomain movesupers(create_movement(tsteps.get_motionstep(), gbxmaps));
  const Observer auto obs = create_observer(config, tsteps, dataset, store);

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

    /* Create zarr store for writing output to storage */
    auto store = FSStore(config.get_zarrbasedir());
    auto dataset = SimpleDataset(store);

    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm = create_sdm(config, tsteps, dataset, store);

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
  std::cout << "-----\n Total Program Duration: " << ttot << "s \n-----\n";

  MPI_Finalize();

  return 0;
}
