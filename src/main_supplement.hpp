// Author: Clara Bayley
// File: main_supplement.hpp
/* File containing functions to create chosen
SDM process and observers to use in main.cpp */

// after make/compiling, execute for example via:
// ./src/coupledCVODECLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"

#ifndef MAIN_SUPPLEMENT_HPP
#define MAIN_SUPPLEMENT_HPP

/* standard library packages */
#include <vector>
#include <stdexcept>
#include <string>

/* constants and equations */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

/* coupled model setup */
#include "./run_coupledmodel_implement.hpp"
#include "./run_coupledmodel.hpp"
#include "./timesteps.hpp"

/* sdm gridboxes setup */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "observers/observers.hpp"
#include "observers/observer_superdropletattributes.hpp"
#include "observers/observer_thermostate.hpp"
#include "observers/zarrstores.hpp"

/* Superdroplet Model (SDM) */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/sdmmotion.hpp"
#include "superdrop_solver/coalescencekernel.hpp"
#include "superdrop_solver/collisionsmethod.hpp"
#include "superdrop_solver/condensationmethod.hpp"
#include "superdrop_solver/sedimentationmethod.hpp"
#include "superdrop_solver/terminalvelocity.hpp"

namespace dlc = dimless_constants;

SdmProcess auto create_sdmprocess(const Config &config,
                                  const ModelTimesteps &mdlsteps)
/* return an SdmProcess type from an amalgamation of other SdmProcess types.
For example return a process that does SDM condensation and collisions from
combined process of those two individual processes */
{
  /* create process for condensation in SDM including Implicit
  Euler Method for solving condensation ODEs */
  const auto condensation_process = CondensationProcess(mdlsteps.condstep,
                                                        &timestep2dimlesstime,
                                                        config.doCouple,
                                                        config.cond_maxiters,
                                                        config.cond_rtol,
                                                        config.cond_atol);

  /* create process for collision-coalescene in SDM */
  const auto probs = GolovinProb(dlc::R0);
  // const auto probs = LongHydrodynamicProb();
  const auto collision_process = CollisionsProcess(mdlsteps.collstep,
                                                   &timestep2realtime,
                                                   probs);

  /* create process for sedimentation in SDM */
  const auto sedimentation_process = SedimentationProcess(mdlsteps.sedistep,
                                                          &timestep2dimlesstime,
                                                          SimmelTerminalVelocity{});
  // const auto sedimentation_process = SedimentationProcess(mdlsteps.sedistep,
  //                                                         RogersYauTerminalVelocity{});

  /* choose an amalgamation of sdm processes to make the returned sdmprocess */
  // const auto sdmprocess =  condensation_process >> collision_process >> sedimentation_process;
  // const auto sdmprocess = condensation_process >> collision_process;
  // const auto sdmprocess = collision_process >> sedimentation_process;
  const auto sdmprocess = collision_process;
  //const auto sdmprocess = condensation_process;
  // const auto sdmprocess = NullProcess{};

  return sdmprocess;
}

SuperdropIntoStoreViaBuffer auto superdropattributes_to_observe()
{
  SuperdropIntoStoreViaBuffer auto id = IdIntoStore();
  SuperdropIntoStoreViaBuffer auto eps = EpsIntoStore();
  SuperdropIntoStoreViaBuffer auto radius = RadiusIntoStore();
  SuperdropIntoStoreViaBuffer auto m_sol = M_solIntoStore();
  SuperdropIntoStoreViaBuffer auto coord3 = Coord3IntoStore();

  SuperdropIntoStoreViaBuffer auto attrs = id >> eps >> radius >> coord3 >> m_sol;

  return attrs;
}

Observer auto create_sdmomentsobserver(SDMomentsStorage &sdmoments)
{
  const Observer auto mom0 = SDMass0thMomentObserver(sdmoments.massmoment0zarr);
  const Observer auto mom1 = SDMassNthMomentObserver(sdmoments.massmoment1zarr, 1);
  const Observer auto mom2 = SDMassNthMomentObserver(sdmoments.massmoment2zarr, 2);

  const auto sdmomentobs = mom2 >> mom1 >> mom0; 

  return sdmomentobs;
}

template <typename ContiguousRaggedZarrStorage>
Observer auto create_observer(const Config &config,
                              ContiguousRaggedZarrStorage &sdzarr,
                              ThermoStateStorage &thermozarr,
                              CoordStorage<double> &timezarr,
                              CoordStorage<unsigned int> &gbxzarr,
                              TwoDStorage<size_t> &nsuperszarr,
                              SDMomentsStorage &sdmoments)
/* return an Observer type from an amalgamation of other observer types.
For example return an observer that observes both the thermostate and the
superdroplets from combination of those two seperate observers */
{
  const Observer auto obs1 = TimeObserver(timezarr);

  const Observer auto obs2 = SDsAttributeObserver(sdzarr);

  const Observer auto obs3 = ThermoStateObserver(thermozarr);
  
  const Observer auto obs4 = GridBoxIndexObserver(gbxzarr);
  
  const Observer auto obs5 = NsupersPerGridBoxObserver(nsuperszarr);

  const Observer auto obs6 = create_sdmomentsobserver(sdmoments);

  const auto observer = obs6 >> obs5 >> obs4 >> obs3 >> obs2 >> obs1 >> PrintObserver{};
  //const auto observer = obs6 >> obs5 >> obs4 >> obs3 >> obs2 >> obs1;

  return observer;
}

#endif // MAIN_SUPPLEMENT_HPP