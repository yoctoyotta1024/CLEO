/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cleo_python_bindings_aliases.hpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Aliases to use to systematically abbreviate various CLEO types,
 * in order to make long template instantiations readable.
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_CLEO_PYTHON_BINDINGS_ALIASES_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_CLEO_PYTHON_BINDINGS_ALIASES_HPP_

#include <pybind11/pybind11.h>

#include "./optional_terminal_velocity.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_motion.hpp"
#include "cartesiandomain/movement/cartesian_transport_across_domain.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "gridboxes/predcorrmotion.hpp"
#include "observers/collect_data_for_simple_dataset.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/gbxindex_observer.hpp"
#include "observers/massmoments_observer.hpp"
#include "observers/nsupers_observer.hpp"
#include "observers/observers.hpp"
#include "observers/sdmmonitor/do_sdmmonitor_obs.hpp"
#include "observers/sdmmonitor/monitor_precipitation_observer.hpp"
#include "observers/state_observer.hpp"
#include "observers/superdrops_observer.hpp"
#include "observers/time_observer.hpp"
#include "observers/totnsupers_observer.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/collisions.hpp"
#include "superdrops/collisions/longhydroprob.hpp"
#include "superdrops/condensation.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "zarr/fsstore.hpp"
#include "zarr/simple_dataset.hpp"

/*
 * aliases as abbreviations of observer types, to make long template of combined observers managable
 */
namespace pyobserver {
using nullmo = NullSDMMonitor;
using precipmo = MonitorPrecipitation;

using gbxindex = GbxindexObserver<SimpleDataset<FSStore>, FSStore>;
using time = ConstTstepObserver<DoTimeObs<SimpleDataset<FSStore>, FSStore>>;
using totnsupers = ConstTstepObserver<DoTotNsupersObs<SimpleDataset<FSStore>, FSStore>>;
using massmoms = ConstTstepObserver<
    DoWriteToDataset<ParallelWriteGridboxes<SimpleDataset<FSStore>, ParallelGridboxesTeamPolicyFunc,
                                            CollectMassMoments<FSStore, MassMomentsFunc>>>>;
using rainmassmoms = ConstTstepObserver<DoWriteToDataset<
    ParallelWriteGridboxes<SimpleDataset<FSStore>, ParallelGridboxesTeamPolicyFunc,
                           CollectMassMoments<FSStore, RaindropsMassMomentsFunc>>>>;
using gridboxes = ConstTstepObserver<DoWriteToDataset<ParallelWriteGridboxes<
    SimpleDataset<FSStore>, ParallelGridboxesRangePolicyFunc,
    CombinedCollectDataForDataset<
        CombinedCollectDataForDataset<
            GenericCollectData<FSStore, uint32_t, NsupersFunc>,
            CombinedCollectDataForDataset<
                CombinedCollectDataForDataset<GenericCollectData<FSStore, float, VvelFunc>,
                                              GenericCollectData<FSStore, float, UvelFunc>>,
                GenericCollectData<FSStore, float, WvelFunc>>>,
        CombinedCollectDataForDataset<
            CombinedCollectDataForDataset<GenericCollectData<FSStore, float, PressFunc>,
                                          GenericCollectData<FSStore, float, TempFunc>>,
            CombinedCollectDataForDataset<GenericCollectData<FSStore, float, QvapFunc>,
                                          GenericCollectData<FSStore, float, QcondFunc>>>>>>>;

using superdrops = ConstTstepObserver<DoWriteToDataset<
    ParallelWriteSupers<SimpleDataset<FSStore>,
                        CombinedCollectDataForDataset<
                            CombinedCollectDataForDataset<
                                CombinedCollectDataForDataset<
                                    CombinedCollectDataForDataset<
                                        CombinedCollectDataForDataset<
                                            CombinedCollectDataForDataset<
                                                CombinedCollectDataForDataset<
                                                    GenericCollectData<FSStore, float, Coord1Func>,
                                                    GenericCollectData<FSStore, float, Coord2Func>>,
                                                GenericCollectData<FSStore, float, Coord3Func>>,
                                            GenericCollectData<FSStore, float, MsolFunc>>,
                                        GenericCollectData<FSStore, float, RadiusFunc>>,
                                    GenericCollectData<FSStore, uint64_t, XiFunc>>,
                                GenericCollectData<FSStore, uint32_t, SdgbxindexFunc>>,
                            GenericCollectData<FSStore, uint32_t, SdIdFunc>>,
                        RaggedCount<SimpleDataset<FSStore>, FSStore>>>>;

using precip = ConstTstepObserver<
    DoSDMMonitorObs<SimpleDataset<FSStore>, FSStore, MonitorPrecipitation, double>>;

using mo01 = CombinedSDMMonitor<nullmo, nullmo>;
using mo012 = CombinedSDMMonitor<mo01, nullmo>;
using mo0123 = CombinedSDMMonitor<mo012, nullmo>;
using mo01234 = CombinedSDMMonitor<mo0123, nullmo>;
using mo012345 = CombinedSDMMonitor<mo01234, nullmo>;
using mo0123456 = CombinedSDMMonitor<mo012345, nullmo>;
using mo01234567 = CombinedSDMMonitor<mo0123456, precipmo>;

using obs01 = CombinedObserver<gbxindex, time, mo01>;
using obs012 = CombinedObserver<obs01, totnsupers, mo012>;
using obs0123 = CombinedObserver<obs012, massmoms, mo0123>;
using obs01234 = CombinedObserver<obs0123, rainmassmoms, mo01234>;
using obs012345 = CombinedObserver<obs01234, gridboxes, mo012345>;
using obs0123456 = CombinedObserver<obs012345, superdrops, mo0123456>;
using obs7 = precip;
using obs = CombinedObserver<obs0123456, obs7, mo01234567>;
}  // namespace pyobserver

/*
 * aliases as abbreviations of types, to make long template instantiations readable.
 * - Abbreviations of concepts/types are as follows:
 *      - map = gridbox maps
 *      - micro = microphysics
 *      - mo = motion
 *      - bcs = boundary conditions
 *      - trans = transport
 *      - move = movement (motion + boundary conditions + transport)
 *      - obs = observer
 * - More specialised abbreviations:
 *      - cart = cartesian
 *      - predcorr = predictor-corrector
 *      - all = SDM with combination of microphysics and superdroplet motion (null observer)
 */
namespace cleo_python_bindings_aliases {
using map_cart = CartesianMaps;

using micro_null = NullMicrophysicalProcess;
using micro_cond = ConstTstepMicrophysics<DoCondensation>;
using micro_colls = ConstTstepMicrophysics<DoCollisions<LongHydroProb, DoCoalescence>>;
using micro_all =
    CombinedMicrophysicalProcess<CombinedMicrophysicalProcess<micro_null, micro_cond>, micro_colls>;

using mo_null = NullMotion;
using mo_cart_predcorr =
    PredCorrMotion<CartesianMaps, OptionalTerminalVelocity, CartesianCheckBounds>;
using bcs_null = NullBoundaryConditions;
using trans_cart = CartesianTransportAcrossDomain;
using move_cart_null = MoveSupersInDomain<map_cart, mo_null, trans_cart, bcs_null>;
using move_cart = MoveSupersInDomain<map_cart, mo_cart_predcorr, trans_cart, bcs_null>;

using obs_null = NullObserver;

using sdm_cart_null = SDMMethods<map_cart, micro_null, mo_null, trans_cart, bcs_null, obs_null>;
using sdm_cart_all =
    SDMMethods<map_cart, micro_all, mo_cart_predcorr, trans_cart, bcs_null, pyobserver::obs>;
}  // namespace cleo_python_bindings_aliases

#endif  // LIBS_CLEO_PYTHON_BINDINGS_CLEO_PYTHON_BINDINGS_ALIASES_HPP_
