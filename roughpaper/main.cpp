/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Monday 29th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 27th March 2024
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

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "coupldyn_fromfile/fromfile_cartesian_dynamics.hpp"
#include "coupldyn_fromfile/fromfilecomms.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "initialise/config.hpp"
#include "initialise/initgbxs_null.hpp"
#include "initialise/initsupers_frombinary.hpp"
#include "initialise/timesteps.hpp"
#include "observers2/observers.hpp"
#include "observers2/streamout_observer.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/initialconditions.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/fsstore.hpp"

template <typename Store>
void test_dataset(Dataset<Store> &dataset) {
  dataset.add_dimension({"SdId", 0});
  auto xzarr = dataset.template create_array<double>("radius", "m", "<f8", 1e-6, {6},
                                                     {"SdId"});  // shape = [0], chunks = 0,1

  dataset.set_dimension({"SdId", 8});
  // dataset.write_to_array(xzarr, h_data);  // shape = [8], chunks = 0,1

  dataset.set_dimension({"SdId", 10});
  dataset.write_arrayshape(xzarr);  // shape = [10], chunks = 0,1
}

/* ---------------------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------------------- */
inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     FSStore &store) {
  const auto obsstep = (unsigned int)tsteps.get_obsstep();
  const auto maxchunk = int{config.maxchunk};

  return StreamOutObserver(obsstep, &step2realtime);
}

inline InitialConditions auto create_initconds(const Config &config) {
  const InitSupersFromBinary initsupers(config);
  const InitGbxsNull initgbxs(config);

  return InitConds(initsupers, initgbxs);
}

inline CoupledDynamics auto create_coupldyn(const Config &config, const CartesianMaps &gbxmaps,
                                            const unsigned int couplstep,
                                            const unsigned int t_end) {
  const auto h_ndims(gbxmaps.ndims_hostcopy());
  const std::array<size_t, 3> ndims({h_ndims(0), h_ndims(1), h_ndims(2)});

  const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);

  return FromFileDynamics(config, couplstep, ndims, nsteps);
}

inline auto create_sdm(const Config &config, const Timesteps &tsteps, FSStore &store) {
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps =
      create_cartesian_maps(config.ngbxs, config.nspacedims, config.grid_filename);
  const MicrophysicalProcess auto microphys = NullMicrophysicalProcess{};
  const Motion<CartesianMaps> auto movesupers = NullMotion{};
  const Observer auto obs = create_observer(config, tsteps, store);
  return SDMMethods(couplstep, gbxmaps, microphys, movesupers, obs);
}

int main(int argc, char *argv[]) {
  Kokkos::Timer kokkostimer;

  /* Read input parameters from configuration file(s) */
  const Config config("/home/m/m300950/CLEO/roughpaper/share/config.txt");
  const Timesteps tsteps(config);  // timesteps for model (e.g. coupling and end time)

  /* Create zarr store for writing output to storage */
  auto fsstore = FSStore(config.zarrbasedir);
  auto dataset = Dataset(fsstore);

  /* Initial conditions for CLEO run */
  const InitialConditions auto initconds = create_initconds(config);

  Kokkos::initialize(argc, argv);
  {
    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm(create_sdm(config, tsteps, fsstore));

    /* Solver of dynamics coupled to CLEO SDM */
    CoupledDynamics auto coupldyn(
        create_coupldyn(config, sdm.gbxmaps, tsteps.get_couplstep(), tsteps.get_t_end()));

    /* coupling between coupldyn and SDM */
    const CouplingComms<FromFileDynamics> auto comms = FromFileComms{};

    /* Run CLEO (SDM coupled to dynamics solver) */
    const RunCLEO runcleo(sdm, coupldyn, comms);
    runcleo(initconds, tsteps.get_t_end());

    test_dataset(dataset);
  }
  Kokkos::finalize();

  const auto ttot = double{kokkostimer.seconds()};
  std::cout << "-------------------------------\n"
               "Total Program Duration: "
            << ttot << "s \n-------------------------------\n";
}
/* ---------------------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------------------- */
