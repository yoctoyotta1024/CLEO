// Author: Clara Bayley
// File: main.cpp
/* This file runs the entire superdrop model (SDM)
thermodynamics read from binaries files
(p, temp, qv and qc) over time for golovin's
collision-coalscence kernel */

// after make/compiling, execute for example via:
// ./src/golcolls0D "../src/config/config.txt" "../libs/claras_SDconstants.hpp"

/* standard library packages */
#include <vector>
#include <stdexcept>
#include <string>

#include <Kokkos_Core.hpp>

/* constants and equations */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

/* sdm gridboxes setup */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/sdmtimesteps.hpp"
#include "sdmgridboxes/runsdmstep.hpp"
#include "sdmgridboxes/sdmotion.hpp"
#include "sdmgridboxes/detectors.hpp"
#include "sdmgridboxes/logbooks.hpp"
#include "observers/observers.hpp"
#include "observers/intostore_observers.hpp"
#include "observers/sdattributes_intostore.hpp"
#include "observers/contigraggedsdstorage.hpp"
#include "observers/thermostatestorage.hpp"
#include "observers/zarrstores.hpp"

/* sdm superdroplets setup */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/coalescencekernel.hpp"
#include "superdrop_solver/collisionsmethod.hpp"

/* thermodynamics solver and coupled model setup */
#include "thermofromfile/run_thermofromfile.hpp"

namespace dlc = dimless_constants;

template <SuperdropIntoStoreViaBuffer S>
struct SomeZarrStores
{
  ThermoStateStorage thermozarr;
  ContiguousRaggedSDStorage<S> sdzarr;
  CoordinateStorage<double> timezarr;

  SomeZarrStores(FSStore &fsstore, const int maxchunk,
                 const unsigned int ngridboxes, S sdattrs)
      : thermozarr(fsstore, maxchunk, ngridboxes),
        sdzarr(fsstore, sdattrs, maxchunk),
        timezarr(fsstore, maxchunk, "time",
                 "<f8", "s", dlc::TIME0) {}
};

SuperdropIntoStoreViaBuffer auto sdattrs_to_observe()
/* choose which methods are used to write attributes
of a superdroplet into zarr storage */
{
  SuperdropIntoStoreViaBuffer auto id = IdIntoStore();
  SuperdropIntoStoreViaBuffer auto eps = EpsIntoStore();
  SuperdropIntoStoreViaBuffer auto radius = RadiusIntoStore();
  SuperdropIntoStoreViaBuffer auto m_sol = M_solIntoStore();

  SuperdropIntoStoreViaBuffer auto attrs = id >> eps >> radius >> m_sol;

  return attrs;
}

template <SuperdropIntoStoreViaBuffer S>
Observer auto create_observer(const int obsstep, SomeZarrStores<S> &stores)
/* return an Observer type from an amalgamation of other observer types.
For example return an observer that observes both the thermostate and the
superdroplets from combination of those two seperate observers */
{
  const Observer auto obs3 = ThermoStateObserver(obsstep, stores.thermozarr);
  const Observer auto obs2 = SDsAttributeObserver(obsstep, stores.sdzarr);
  const Observer auto obs1 = TimeObserver(obsstep, stores.timezarr);

  const auto observer = obs3 >> obs2 >> obs1 >> PrintObserver(obsstep);

  return observer;
}

int main(int argc, char *argv[])
{
  Kokkos::Timer kokkostimer;

  if (argc < 3)
  {
    throw std::invalid_argument("config and/or constants files not specified");
  }

  /* object containing input parameters from configuration file */
  const std::string configfilepath = argv[1];    // path to configuration (.txt file)
  const std::string constantsfilepath = argv[2]; // path to constants (.hpp file)
  const Config config(configfilepath, constantsfilepath);

  /* object for time-stepping parameters of coupled model */
  const SDMTimesteps mdlsteps(config.CONDTSTEP, config.COLLTSTEP,
                              config.MOTIONTSTEP, config.COUPLTSTEP,
                              config.OBSTSTEP, config.T_END);

  /* create map from gridbox index to its coordinate boundaries */
  const Maps4GridBoxes gbxmaps(config.SDnspace, config.grid_filename);

  /* create superdroplet model (SDM) process from combination of chosen SDM processes */
  const auto sdmprocess(CollisionsProcess(mdlsteps.collsubstep,
                                          &step2realtime,
                                          GolovinProb(dlc::R0)));
  const MoveSuperdropsInDomain sdmmotion(NullMotion{});

  /* create observer from combination of chosen observers */
  FSStore fsstore(config.zarrbasedir);
  SomeZarrStores zarrstores(fsstore, config.maxchunk,
                            gbxmaps.gbxidxs.size(),
                            sdattrs_to_observe());
  const auto observer = create_observer(mdlsteps.obsstep, zarrstores);

  // DetectionLogbooks logbooks;
  // const DetectorsInstallation dtrs(logbooks);

  const RunSDMStep sdm(gbxmaps, sdmmotion, sdmprocess, observer);

  Kokkos::initialize(argc, argv);
  {
    /* RUN SDM MODEL WITH THERMODYNAMICS FROM FILE */
    run_thermofromfile(config, sdm, dtrs, mdlsteps.t_end, mdlsteps.couplstep);
  }
  Kokkos::finalize();
  std::cout << "  ------ Total Duration: " << kokkostimer.seconds() << "s ----- \n";

  return 0;
}