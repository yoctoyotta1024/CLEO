/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_impl.hpp
 * Project: src
 * Created Date: Monday 29th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for main.cpp to run the CLEO super-droplet model (SDM).
 */

#ifndef ROUGHPAPER_SRC_MAIN_IMPL_HPP_
#define ROUGHPAPER_SRC_MAIN_IMPL_HPP_

#include <Kokkos_Core.hpp>
#include <cmath>
#include <concepts>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string_view>

#include "cartesiandomain/add_supers_at_domain_top.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/cartesianmotion.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"
#include "cartesiandomain/null_boundary_conditions.hpp"
#include "coupldyn_fromfile/fromfile_cartesian_dynamics.hpp"
#include "coupldyn_fromfile/fromfilecomms.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "initialise/config.hpp"
#include "initialise/init_all_supers_from_binary.hpp"
#include "initialise/init_supers_from_binary.hpp"
#include "initialise/initgbxsnull.hpp"
#include "initialise/initialconditions.hpp"
#include "initialise/timesteps.hpp"
#include "observers/gbxindex_observer.hpp"
#include "observers/massmoments_observer.hpp"
#include "observers/nsupers_observer.hpp"
#include "observers/observers.hpp"
#include "observers/runstats_observer.hpp"
#include "observers/sdmmonitor/monitor_condensation_observer.hpp"
#include "observers/sdmmonitor/monitor_massmoments_observer.hpp"
#include "observers/state_observer.hpp"
#include "observers/streamout_observer.hpp"
#include "observers/superdrops_observer.hpp"
#include "observers/time_observer.hpp"
#include "observers/totnsupers_observer.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/collisions/breakup.hpp"
#include "superdrops/collisions/breakup_nfrags.hpp"
#include "superdrops/collisions/coalbure.hpp"
#include "superdrops/collisions/coalbure_flag.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/constprob.hpp"
#include "superdrops/collisions/golovinprob.hpp"
#include "superdrops/collisions/longhydroprob.hpp"
#include "superdrops/collisions/lowlistprob.hpp"
#include "superdrops/condensation.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/terminalvelocity.hpp"
#include "zarr/dataset.hpp"
#include "zarr/fsstore.hpp"

inline CoupledDynamics auto create_coupldyn(const Config &config, const CartesianMaps &gbxmaps,
                                            const unsigned int couplstep,
                                            const unsigned int t_end) {
  const auto h_ndims(gbxmaps.ndims_hostcopy());
  const std::array<size_t, 3> ndims({h_ndims(0), h_ndims(1), h_ndims(2)});

  const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);

  return FromFileDynamics(config.get_fromfiledynamics(), couplstep, ndims, nsteps);
}

inline InitialConditions auto create_initconds(const Config &config) {
  // const InitAllSupersFromBinary initsupers(config.get_initsupersfrombinary());
  const InitSupersFromBinary initsupers(config.get_initsupersfrombinary());
  const InitGbxsNull initgbxs(config.get_ngbxs());

  return InitConds(initsupers, initgbxs);
}

inline GridboxMaps auto create_gbxmaps(const Config &config) {
  const auto gbxmaps = create_cartesian_maps(config.get_ngbxs(), config.get_nspacedims(),
                                             config.get_grid_filename());
  return gbxmaps;
}

inline Motion<CartesianMaps> auto create_motion(const unsigned int motionstep) {
  // const auto terminalv = NullTerminalVelocity{};
  // const auto terminalv = RogersYauTerminalVelocity{};
  // const auto terminalv = SimmelTerminalVelocity{};
  const auto terminalv = RogersGKTerminalVelocity{};

  return CartesianMotion(motionstep, &step2dimlesstime, terminalv);

  // return NullMotion{};
}

inline auto create_boundary_conditions(const Config &config) {
  // return AddSupersAtDomainTop(config.get_addsupersatdomaintop());
  return NullBoundaryConditions{};
}

template <GridboxMaps GbxMaps>
inline auto create_movement(const Config &config, const Timesteps &tsteps, const GbxMaps &gbxmaps) {
  const auto motion(create_motion(tsteps.get_motionstep()));

  const auto boundary_conditions(create_boundary_conditions(config));

  return MoveSupersInDomain(gbxmaps, motion, boundary_conditions);
}

inline MicrophysicalProcess auto config_condensation(const Config &config,
                                                     const Timesteps &tsteps) {
  const auto c = config.get_condensation();

  return Condensation(tsteps.get_condstep(), &step2dimlesstime, c.do_alter_thermo, c.maxniters,
                      c.rtol, c.atol, c.MINSUBTSTEP, &realtime2dimless);
}

inline MicrophysicalProcess auto config_collisions(const Config &config, const Timesteps &tsteps) {
  // const PairProbability auto collprob = LongHydroProb();
  // // const NFragments auto nfrags = ConstNFrags(config.get_breakup().constnfrags.nfrags);
  // const NFragments auto nfrags = CollisionKineticEnergyNFrags{};
  // // const CoalBuReFlag auto coalbure_flag = SUCoalBuReFlag{};
  // const CoalBuReFlag auto coalbure_flag = TSCoalBuReFlag{};
  // const MicrophysicalProcess auto colls = CoalBuRe(tsteps.get_collstep(),
  //                                                  &step2realtime,
  //                                                  collprob,
  //                                                  nfrags,
  //                                                  coalbure_flag);
  // return colls;

  // const PairProbability auto buprob = LowListBuProb();
  // const NFragments auto nfrags = ConstNFrags(config.get_breakup().constnfrags.nfrags);
  // const MicrophysicalProcess auto bu = CollBu(tsteps.get_collstep(),
  //                                             &step2realtime,
  //                                             buprob,
  //                                             nfrags);

  // const PairProbability auto coalprob = LowListCoalProb();
  // const PairProbability auto coalprob = GolovinProb();
  const PairProbability auto coalprob = LongHydroProb(1.0);
  const MicrophysicalProcess auto coal = CollCoal(tsteps.get_collstep(), &step2realtime, coalprob);

  return coal;
  // return coal >> bu;
}

inline MicrophysicalProcess auto create_microphysics(const Config &config,
                                                     const Timesteps &tsteps) {
  const MicrophysicalProcess auto cond = config_condensation(config, tsteps);
  // const MicrophysicalProcess auto colls = config_collisions(config, tsteps);
  // return colls >> cond;

  // const MicrophysicalProcess auto null = NullMicrophysicalProcess{};
  // return null;

  return cond;
}

template <typename Store>
inline Observer auto create_superdrops_observer(const unsigned int interval,
                                                Dataset<Store> &dataset, const size_t maxchunk) {
  CollectDataForDataset<Store> auto sdid = CollectSdId(dataset, maxchunk);
  CollectDataForDataset<Store> auto sdgbxindex = CollectSdgbxindex(dataset, maxchunk);
  CollectDataForDataset<Store> auto xi = CollectXi(dataset, maxchunk);
  CollectDataForDataset<Store> auto radius = CollectRadius(dataset, maxchunk);
  CollectDataForDataset<Store> auto msol = CollectMsol(dataset, maxchunk);
  CollectDataForDataset<Store> auto coord3 = CollectCoord3(dataset, maxchunk);
  // CollectDataForDataset<Store> auto coord1 = CollectCoord1(dataset, maxchunk);
  // CollectDataForDataset<Store> auto coord2 = CollectCoord2(dataset, maxchunk);

  const auto collect_sddata = coord3 >> msol >> radius >> xi >> sdgbxindex >> sdid;
  return SuperdropsObserver(interval, dataset, maxchunk, collect_sddata);
}

template <typename Store>
inline Observer auto create_gridboxes_observer(const unsigned int interval, Dataset<Store> &dataset,
                                               const size_t maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Store> auto thermo = CollectThermo(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Store> auto windvel = CollectWindVel(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Store> auto nsupers = CollectNsupers(dataset, maxchunk, ngbxs);

  const CollectDataForDataset<Store> auto collect_gbxdata = nsupers >> windvel >> thermo;
  return WriteToDatasetObserver(interval, dataset, collect_gbxdata);
}

template <typename Store>
inline Observer auto create_sdmmonitor_observer(const unsigned int interval,
                                                Dataset<Store> &dataset, const size_t maxchunk,
                                                const size_t ngbxs) {
  const Observer auto obs_cond = MonitorCondensationObserver(interval, dataset, maxchunk, ngbxs);
  const Observer auto obs_massmoms = MonitorMassMomentsObserver(interval, dataset, maxchunk, ngbxs);

  return obs_cond >> obs_massmoms;
}

template <typename Store>
inline Observer auto create_observer(const Config &config, const Timesteps &tsteps,
                                     Dataset<Store> &dataset) {
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();
  const auto ngbxs = config.get_ngbxs();

  const Observer auto obs0 = RunStatsObserver(obsstep, config.get_stats_filename());

  const Observer auto obs1 = StreamOutObserver(obsstep, &step2realtime);

  const Observer auto obs2 = TimeObserver(obsstep, dataset, maxchunk, &step2dimlesstime);

  const Observer auto obs3 = GbxindexObserver(dataset, maxchunk, ngbxs);

  const Observer auto obs4 = TotNsupersObserver(obsstep, dataset, maxchunk);

  const Observer auto obs5 = MassMomentsObserver(obsstep, dataset, maxchunk, ngbxs);

  const Observer auto obs6 = MassMomentsRaindropsObserver(obsstep, dataset, maxchunk, ngbxs);

  const Observer auto obsgbx = create_gridboxes_observer(obsstep, dataset, maxchunk, ngbxs);

  const Observer auto obssd = create_superdrops_observer(obsstep, dataset, maxchunk);

  const Observer auto obsm = create_sdmmonitor_observer(obsstep, dataset, maxchunk, ngbxs);

  return obsm >> obssd >> obsgbx >> obs6 >> obs5 >> obs4 >> obs3 >> obs2 >> obs1 >> obs0;
}

template <typename Store>
inline auto create_sdm(const Config &config, const Timesteps &tsteps, Dataset<Store> &dataset) {
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps(create_gbxmaps(config));
  const MicrophysicalProcess auto microphys(create_microphysics(config, tsteps));
  const MoveSupersInDomain movesupers(create_movement(config, tsteps, gbxmaps));
  const Observer auto obs(create_observer(config, tsteps, dataset));

  return SDMMethods(couplstep, gbxmaps, microphys, movesupers, obs);
}

#endif  // ROUGHPAPER_SRC_MAIN_IMPL_HPP_
