/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_fromfile.cpp
 * Project: src
 * Created Date: Friday 22nd March 2024
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
 * runs the CLEO super-droplet model (SDM) for 3-D setup reading data from file
 * after make/compiling, execute for example via:
 * ./src/fromfile ../src/config/config.yaml
 */

#include <Kokkos_Core.hpp>
#include <cmath>
#include <concepts>
#include <iostream>
#include <stdexcept>
#include <string_view>

#include "cartesiandomain/cartesian_decomposition.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_motion.hpp"
#include "cartesiandomain/movement/cartesian_movement.hpp"
#include "configuration/config.hpp"
#include "coupldyn_fromfile/fromfile_cartesian_dynamics.hpp"
#include "coupldyn_fromfile/fromfilecomms.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "initialise/init_supers_from_binary.hpp"
#include "initialise/initgbxsnull.hpp"
#include "initialise/initialconditions.hpp"
#include "initialise/timesteps.hpp"
#include "mpi.h"
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
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "zarr/collective_dataset.hpp"
#include "zarr/fsstore.hpp"

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
  const auto initsupers = InitSupersFromBinary(config.get_initsupersfrombinary(), gbxmaps);
  const auto initgbxs = InitGbxsNull(gbxmaps.get_local_ngridboxes_hostcopy());

  return InitConds(initsupers, initgbxs);
}

inline GridboxMaps auto create_gbxmaps(const Config &config) {
  const auto gbxmaps = create_cartesian_maps(config.get_ngbxs(), config.get_nspacedims(),
                                             config.get_grid_filename());
  return gbxmaps;
}

inline MicrophysicalProcess auto create_microphysics(const Config &config,
                                                     const Timesteps &tsteps) {
  return NullMicrophysicalProcess{};
}

inline auto create_movement(const unsigned int motionstep, const CartesianMaps &gbxmaps) {
  const auto terminalv = NullTerminalVelocity{};
  const Motion<CartesianMaps> auto motion =
      CartesianMotion(motionstep, &step2dimlesstime, terminalv);

  const BoundaryConditions<CartesianMaps> auto boundary_conditions = NullBoundaryConditions{};

  return cartesian_movement(gbxmaps, motion, boundary_conditions);
}

template <typename Store>
inline Observer auto create_superdrops_observer(const unsigned int interval,
                                                Dataset<Store> &dataset, const int maxchunk) {
  CollectDataForDataset<Store> auto sdid = CollectSdId(dataset, maxchunk);
  CollectDataForDataset<Store> auto coord3 = CollectCoord3(dataset, maxchunk);
  CollectDataForDataset<Store> auto coord1 = CollectCoord1(dataset, maxchunk);
  CollectDataForDataset<Store> auto coord2 = CollectCoord2(dataset, maxchunk);

  const auto collect_data = sdid >> coord3 >> coord1 >> coord2;
  return SuperdropsObserver(interval, dataset, maxchunk, collect_data);
}

template <typename Store>
inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     Dataset<Store> &dataset, const CartesianMaps &gbxmaps) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();
  const auto ngbxs = gbxmaps.get_local_ngridboxes();

  const Observer auto obs0 = StreamOutObserver(obsstep, &step2realtime);

  const Observer auto obs1 = TimeObserver(obsstep, dataset, maxchunk, &step2dimlesstime);

  const Observer auto obs2 = GbxindexObserver(dataset, maxchunk, ngbxs);

  const Observer auto obs3 =
      StateObserver(obsstep, dataset, maxchunk, gbxmaps.get_total_global_ngridboxes());

  const Observer auto obssd = create_superdrops_observer(obsstep, dataset, maxchunk);

  return obssd >> obs3 >> obs2 >> obs1 >> obs0;
}

template <typename Store>
inline auto create_sdm(const Config &config, const Timesteps &tsteps, Dataset<Store> &dataset) {
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps(create_gbxmaps(config));
  const MicrophysicalProcess auto microphys(create_microphysics(config, tsteps));
  const MoveSupersInDomain movesupers(create_movement(tsteps.get_motionstep(), gbxmaps));
  const Observer auto obs(create_observer(config, tsteps, dataset, gbxmaps));

  return SDMMethods(couplstep, gbxmaps, microphys, movesupers, obs);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    throw std::invalid_argument("configuration file(s) not specified");
  }

  MPI_Init(&argc, &argv);

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

    /* Adjust dataset given domain decomposition */
    dataset.set_decomposition(sdm.gbxmaps.get_domain_decomposition());
    dataset.set_max_superdroplets(config.get_maxnsupers());

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
