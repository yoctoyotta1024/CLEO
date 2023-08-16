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
#include "sdmgridboxes/detectors.hpp"
#include "sdmgridboxes/detectors_ptr.hpp"
#include "sdmgridboxes/logbooks.hpp"

/* sdm observers setup */
#include "observers/observers.hpp"
#include "observers/observegbxs.hpp"
#include "observers/observelbks.hpp"
#include "observers/gridboxes_intostore.hpp"
#include "observers/logbooks_intostore.hpp"

#include "singlevarstorage.hpp"
#include "massmomentsstorage.hpp"
#include "zarrstorage/sdattributes_intostore.hpp"
#include "zarrstorage/contigraggedsdstorage.hpp"
#include "zarrstorage/thermostatestorage.hpp"
#include "zarrstorage/zarrstores.hpp"

/* sdm superdroplets setup */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/coalescence.hpp"
#include "superdrop_solver/breakup.hpp"
#include "superdrop_solver/coal_breakup_rebound.hpp"
#include "superdrop_solver/coal_breakup.hpp"
#include "superdrop_solver/condensation.hpp"
#include "superdrop_solver/sedimentation.hpp"
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
  MomentsStorages massmomszarr;
  RainMomentsStorages rainmassmomszarr;
  CoordinateStorage<double> timezarr;
  CoordinateStorage<unsigned int> gbxzarr;
  TwoDStorage<size_t> nsuperszarr;
  TwoDStorage<size_t> nrainsuperszarr;
  LogbooksStorage<double> lbkszarr;

  SomeZarrStores(FSStore &store, const int maxchunk,
                 const unsigned int ngbxs, S sdattrs)
      : thermozarr(store, maxchunk, ngbxs),
        sdzarr(store, sdattrs, maxchunk),
        sdgbxzarr(store, SdgbxIntoStore(), maxchunk),
        massmomszarr(store, maxchunk, ngbxs),
        rainmassmomszarr(store, maxchunk, ngbxs),
        timezarr(make_timezarr(store, maxchunk)),
        gbxzarr(make_gbxzarr(store, maxchunk)),
        nsuperszarr(make_nsuperszarr(store, maxchunk, ngbxs)),
        nrainsuperszarr(make_nrainsuperszarr(store, maxchunk, ngbxs)),
        lbkszarr(make_logbookszarr(store, maxchunk)) {}
};

SdMotion auto create_sdmotion(const int motionstep)
{
  // const auto terminalv = RogersYauTerminalVelocity{};
  // const auto terminalv = NullTerminalVelocity{};
  const auto terminalv = SimmelTerminalVelocity{};
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
  const double cond_subtstep(realtime2dimless(config.cond_SUBTSTEP));
  const auto cond(CondensationProcess(mdlsteps.condsubstep,
                                      &step2dimlesstime,
                                      config.doAlterThermo,
                                      config.cond_iters,
                                      cond_subtstep,
                                      config.cond_rtol,
                                      config.cond_atol));

  /* create process for collision-coalescene in SDM */
  // const auto terminalv(SimmelTerminalVelocity{});
  // const auto probs_coal(CollCoalProb_LowList(terminalv));
  // const auto probs_coal = CollCoalProb_Golovin();
  const auto probs_coal(CollCoalProb_Long());
  const auto coal(CollisionCoalescenceProcess(mdlsteps.collsubstep,
                                              &step2realtime,
                                              probs_coal));

  // /* create process for collision-breakup in SDM */
  // const auto terminalv = SimmelTerminalVelocity{};
  // const auto probs_bu = CollBuProb_LowList(terminalv);
  // const auto bu(CollisionBreakupProcess(mdlsteps.collsubstep,
  //                                       &step2realtime,
  //                                       probs_bu,
  //                                       config.nfrags));

  // const auto probs_coll(CollCoalProb_Long());
  // const auto terminalv(SimmelTerminalVelocity{});
  // const auto collall(CollisionAllProcess(mdlsteps.collsubstep,
  //                                        &step2realtime,
  //                                        probs_coll,
  //                                        terminalv,
  //                                        config.nfrags));

  // const double coalrate(5e-7); // coalescence rate [s^-1]
  // const double burate(0.0);  // breakup rate [s^-1]
  // const auto collcnst(CollisionCoalBuConst(mdlsteps.collsubstep,
  //                                         &step2realtime,
  //                                         config.nfrags,
  //                                         coalrate, burate));

  /* choose an amalgamation of sdm processes to make the returned sdmprocess */
  const auto sdmprocess = cond >> coal;
  // const auto sdmprocess = coal;
  // const auto sdmprocess = cond;
  // const auto sdmprocess = coalall;
  // const auto sdmprocess = collcnst;

  // return NullProcess{};
  return sdmprocess;
}

CreateDetectorsPtr auto specify_detectors(const Maps4GridBoxes &gbxmaps)
{
  // return NullDetectorsPtr{};
  return PrecipDetectorsPtr(gbxmaps);
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

  SuperdropIntoStoreViaBuffer auto attrs = id >> eps >> radius >>
                                           m_sol >> coord3 >> coord1;

  return attrs;
}

ObserveGBxs auto
create_observegbx_massmoments(MomentsStorages &mms,
                              RainMomentsStorages &rmms,
                              const size_t ngbxs)
{
  const auto mom0 = ObserveNthMassMoment(mms.mom0zarr, 0, ngbxs);
  const auto mom1 = ObserveNthMassMoment(mms.mom1zarr, 1, ngbxs);
  const auto mom2 = ObserveNthMassMoment(mms.mom2zarr, 2, ngbxs);

  const auto rain0 = ObserveNthRainMassMoment(rmms.mom0zarr, 0, ngbxs);
  const auto rain1 = ObserveNthRainMassMoment(rmms.mom1zarr, 1, ngbxs);
  // const auto rain2 = ObserveNthRainMassMoment(rmms.mom2zarr, 2, ngbxs);
  
  return rain1 >> rain0 >> mom2 >> mom1 >> mom0; 
}

template <SuperdropIntoStoreViaBuffer S>
Observer auto create_observer(SomeZarrStores<S> &stores, const int obsstep,
                              const size_t ngbxs)
/* return an Observer type from an amalgamation of other observer types.
For example return an observer that observes both the thermostate and the
superdroplets from combination of those two seperate observers */
{
  const ObserveGBxs auto og1 = ObserveTime(stores.timezarr);

  const ObserveGBxs auto og2a = ObserveSDsAttributes(stores.sdzarr);
  const ObserveGBxs auto og2b = ObserveSDsGbxindex(stores.sdgbxzarr);

  // const ObserveGBxs auto og3 = ObserveThermoState(stores.thermozarr);
  
  const ObserveGBxs auto og4 = ObserveGridBoxIndex(stores.gbxzarr);

  const ObserveGBxs auto og5 = ObserveNsupersPerGridBox(stores.nsuperszarr,
                                                        ngbxs);

  const ObserveGBxs auto og6 = ObserveNRainsupersPerGridBox(stores.nrainsuperszarr,
                                                            ngbxs);

  const ObserveGBxs auto og7 = create_observegbx_massmoments(stores.massmomszarr,
                                                            stores.rainmassmomszarr,
                                                            ngbxs);

  const ObserveLbks auto ol1 = ObservePrecip(stores.lbkszarr);

  // const ObserveGBxs auto obsgbxs = og6 >> og5 >> og4 >> og3 >> og2a >> og2b >> og1;
  const ObserveGBxs auto obsgbxs = og7 >> og6 >> og5 >> og4 >> og2a >> og2b >> og1;

  const Observer auto obs1 = ConstIntervalGBxsObserver(obsstep, obsgbxs);
  const Observer auto obs2 = ConstIntervalLbksObserver(obsstep, ol1);

  const Observer auto observer = obs1 >> obs2 >> PrintObserver(obsstep);
  // const Observer auto observer = obs1 >> obs2;

  return observer;
}

#endif // MAIN_SUPPLEMENT_HPP