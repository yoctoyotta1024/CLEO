// Author: Clara Bayley
// File: main_supplement.hpp
/* File containing functions to create chosen
SDM process and observers to use in main.cpp */

// after make/compiling, execute for example via:
// ./src/runCLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"

#ifndef MAIN_SUPPLEMENT_HPP
#define MAIN_SUPPLEMENT_HPP

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
#include "superdrop_solver/condensationmethod.hpp"
#include "superdrop_solver/sedimentationmethod.hpp"
#include "superdrop_solver/terminalvelocity.hpp"

/* thermodynamics solver and coupled model setup */
#include "thermofromfile/run_thermofromfile.hpp"
#include "thermofromfile/prescribedmotion.hpp"

namespace dlc = dimless_constants;

template <SuperdropIntoStoreViaBuffer S>
struct SomeZarrStores
{
  ThermoStateStorage thermozarr;
  ContiguousRaggedSDStorage<S> sdzarr;
  ContiguousRaggedSDStorage<SdgbxIntoStore> sdgbxzarr;
  MassMomStorages massmoms;
  RainMassMomStorages rainmassmoms;
  CoordinateStorage<double> timezarr;
  CoordinateStorage<unsigned int> gbxzarr;
  TwoDStorage<size_t> nsuperszarr;

SomeZarrStores(FSStore &fsstore, const int maxchunk,
              const unsigned int ngridboxes, S sdattrs)
      : thermozarr(fsstore, maxchunk, ngridboxes),
        sdzarr(fsstore, sdattrs, maxchunk),
        sdgbxzarr(fsstore, SdgbxIntoStore(), maxchunk),
        massmoms(fsstore, maxchunk, ngridboxes),
        rainmassmoms(fsstore, maxchunk, ngridboxes),
        timezarr(fsstore, maxchunk, "time",
                 "<f8", "s", dlc::TIME0),
        gbxzarr(fsstore, maxchunk, "gbxindex",
                "<u4", " ", 1),
        nsuperszarr(fsstore, maxchunk, "nsupers",
                    "<u8", " ", 1, ngridboxes) {}
};

SdMotion auto create_sdmotion(const int motionstep)
{
  // const auto terminalv = RogersYauTerminalVelocity{};
  const auto terminalv = SimmelTerminalVelocity{};
  // const auto terminalv = NullTerminalVelocity{};
  const SdMotion auto movewithsedi = MoveWithSedimentation(motionstep,
                                                          &step2dimlesstime,
                                                          terminalv);

  // auto rhotilda = [](const ThermoState &state)
  // { return state.press / (state.temp * (dlc::Rgas_dry + state.qvap * dlc::Rgas_v)); };
  // const Prescribed2DFlow flow2d(1500 / dlc::COORD0, 1500 / dlc::COORD0,
  //                               0.6 / dlc::W0, rhotilda);
  // const SdMotion auto prescribed2d = MoveWith2DPrescribedFlow(motionstep,
  //                                                 &step2dimlesstime,
  //                                                 flow2d);
  
  return movewithsedi;
  // return prescribed2d;
  // return NullMotion{};
}

SdmProcess auto create_sdmprocess(const Config &config,
                                  const SDMTimesteps &mdlsteps)
/* return an SdmProcess type from an amalgamation of other SdmProcess types.
For example return a process that does SDM condensation and collisions from
combined process of those two individual processes */
{
  /* create process for condensation in SDM including Implicit
  Euler Method for solving condensation ODEs */
  const double cond_subtstep = realtime2dimless(config.cond_SUBTSTEP);
  const auto cond(CondensationProcess(mdlsteps.condsubstep, &step2dimlesstime,
                                            config.doAlterThermo, config.cond_iters,
                                            cond_subtstep, config.cond_rtol,
                                            config.cond_atol));

  /* create process for collision-coalescene in SDM */
  const auto probs(GolovinProb(dlc::R0));
  // const auto probs(LongHydrodynamicProb());
  const auto colls(CollisionsProcess(mdlsteps.collsubstep,
                                       &step2realtime,
                                       probs));

  /* create process for sedimentation in SDM -> n.b. this has moved to sdmmotion*/
  // // const auto terminalv(RogersYauTerminalVelocity{});
  // const auto terminalv(SimmelTerminalVelocity{});
  // const auto sedi(SedimentationProcess(mdlsteps.motionstep,
  //                                        &step2dimlesstime,
  //                                        terminalv));

  /* choose an amalgamation of sdm processes to make the returned sdmprocess */
  const auto sdmprocess = cond >> colls;
  // const auto sdmprocess = cond;
  // const auto sdmprocess = colls;

  return sdmprocess;
  // return NullProcess{};
}

SuperdropIntoStoreViaBuffer auto sdattrs_to_observe()
/* choose which methods are used to write attributes
of a superdroplet into zarr storage */
{
  SuperdropIntoStoreViaBuffer auto id = IdIntoStore();
  SuperdropIntoStoreViaBuffer auto eps = EpsIntoStore();
  SuperdropIntoStoreViaBuffer auto radius = RadiusIntoStore();
  SuperdropIntoStoreViaBuffer auto m_sol = M_solIntoStore();
  SuperdropIntoStoreViaBuffer auto coord3 = Coord3IntoStore();
  SuperdropIntoStoreViaBuffer auto coord1 = Coord1IntoStore();
  // SuperdropIntoStoreViaBuffer auto coord2 = Coord2IntoStore();

  SuperdropIntoStoreViaBuffer auto attrs = id >> eps >> radius >> m_sol >> coord3 >> coord1;

  return attrs;
}

Observer auto create_massmoments_observer(const int obsstep,
                                          MassMomStorages &mms,
                                          RainMassMomStorages &rmms)
{
  const Observer auto mom0 = NthMassMomentObserver(obsstep, mms.mom0zarr, 0);
  const Observer auto mom1 = NthMassMomentObserver(obsstep, mms.mom1zarr, 1);
  const Observer auto mom2 = NthMassMomentObserver(obsstep, mms.mom2zarr, 2);

  const Observer auto rain0 = NthRainMassMomentObserver(obsstep, rmms.mom0zarr, 0);
  const Observer auto rain1 = NthRainMassMomentObserver(obsstep, rmms.mom1zarr, 1);
  // const Observer auto rain2 = NthRainMassMomentObserver(obsstep, rmms.mom2zarr, 2);
  
  const auto momsobs = rain1 >> rain0 >> mom2 >> mom1 >> mom0; 

  return momsobs;
}

template <SuperdropIntoStoreViaBuffer S>
Observer auto create_observer(const int obsstep, SomeZarrStores<S> &stores)
/* return an Observer type from an amalgamation of other observer types.
For example return an observer that observes both the thermostate and the
superdroplets from combination of those two seperate observers */
{
  const Observer auto obs1 = TimeObserver(obsstep, stores.timezarr);

  const Observer auto obs2a = SDsAttributeObserver(obsstep, stores.sdzarr);
  const Observer auto obs2b = SDsGbxindexObserver(obsstep, stores.sdgbxzarr);

  // const Observer auto obs3 = ThermoStateObserver(obsstep, stores.thermozarr);
  
  const Observer auto obs4 = GridBoxIndexObserver(obsstep, stores.gbxzarr);
  
  const Observer auto obs5 = NsupersPerGridBoxObserver(obsstep, stores.nsuperszarr);

  const Observer auto obs6 = create_massmoments_observer(obsstep,
                                                         stores.massmoms,
                                                         stores.rainmassmoms);

  // const auto observer = obs6 >> obs5 >> obs4 >> obs3 >>
  //                       obs2a >> obs2b >> obs1 >> PrintObserver(obsstep);
  const auto observer = obs6 >> obs5 >> obs4 >>
                        obs2a >> obs2b >> obs1;

  return observer;
}

#endif // MAIN_SUPPLEMENT_HPP