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

/* constants and equations */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

/* coupled model setup */
#include "./run_coupledmodel_implement.hpp"
#include "./run_coupledmodel.hpp"
#include "./runmodel/maps4gridboxes.hpp"
#include "./runmodel/timesteps.hpp"
#include "./runmodel/observers/observers.hpp"
#include "./runmodel/observers/observer_superdropletattributes.hpp"
#include "./runmodel/observers/observer_thermostate.hpp"
#include "./runmodel/observers/zarrstores.hpp"

/* Superdroplet Model (SDM) */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/coalescencekernel.hpp"
#include "superdrop_solver/collisionsmethod.hpp"
#include "superdrop_solver/condensationmethod.hpp"
#include "superdrop_solver/sedimentationmethod.hpp"
#include "superdrop_solver/terminalvelocity.hpp"

namespace dlc = dimless_constants;

SdmProcess auto create_sdmprocess(const Config &config,
                                  const Timesteps &mdlsteps)
/* return an SdmProcess type from an amalgamation of other SdmProcess types.
For example return a process that does SDM condensation and collisions from
combined process of those two individual processes */
{
  /* create process for collision-coalescene in SDM */
  const auto probs = GolovinProb(dlc::R0);
  //const auto probs = LongHydrodynamicProb();
  const double COLLTSTEP = timestep2realtime(mdlsteps.collstep);
  const CollisionsMethod collisions(COLLTSTEP, probs);
  const auto collision_process = ConstTstepProcess{mdlsteps.collstep,
                                                   collisions};

  /* create process for condensation in SDM including Implicit
  Euler Method for solving condensation ODEs */
  const double condtimestep = timestep2dimlesstime(mdlsteps.condstep);
  const CondensationMethod condensationmethod(config.doCouple, condtimestep,
                                              config.cond_maxiters, config.cond_rtol,
                                              config.cond_atol);
  const auto condensation_process = ConstTstepProcess{mdlsteps.condstep,
                                                      condensationmethod};

  /* create process for sedimentation in SDM */
  const double seditstep = timestep2dimlesstime(mdlsteps.sedistep);
  const SedimentationMethod sedimentation(seditstep, SimmelTerminalVelocity{});
  // const SedimentationMethod sedimentation(seditstep, RogersYauTerminalVelocity{});
  const auto sedimentation_process = ConstTstepProcess{mdlsteps.sedistep,
                                                       sedimentation};

  /* choose an amalgamation of sdm processes to make the returned sdmprocess */
  // const auto sdmprocess =  condensation_process >> collision_process >> sedimentation_process;
  // const auto sdmprocess = condensation_process >> collision_process;
  // const auto sdmprocess = collision_process >> sedimentation_process;
  const auto sdmprocess = collision_process;
  // const auto sdmprocess = condensation_process;
  // const auto sdmprocess = NullProcess{};

  return sdmprocess;
}

template <typename ContiguousRaggedZarrStorage>
Observer auto create_observer(const Config &config,
                              ContiguousRaggedZarrStorage &sdzarr,
                              ThermoStateStorage &thermozarr,
                              CoordStorage<double> &timezarr,
                              CoordStorage<unsigned int> &gbxzarr,
                              TwoDStorage<size_t> &nsuperszarr)
/* return an Observer type from an amalgamation of other observer types.
For example return an observer that observes both the thermostate and the
superdroplets from combination of those two seperate observers */
{
  const Observer auto obs1 = TimeObserver(timezarr);

  const Observer auto obs2 = SDsAttributeObserver(sdzarr);

  const Observer auto obs3 = ThermoStateObserver(thermozarr);
  
  const Observer auto obs4 = GridBoxIndexObserver(gbxzarr);
  
  const Observer auto obs5 = NsupersPerGridBoxObserver(nsuperszarr);

  const auto observer = obs5 >> obs4 >> obs3 >> obs2 >> obs1 >> PrintObserver{};

  return observer;
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

#endif // MAIN_SUPPLEMENT_HPP