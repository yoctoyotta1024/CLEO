/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: pycleo_aliases.hpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 15th July 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Aliases to use to systematically abbreviate various CLEO types,
 * in order to make long template instantiations readable.
 */

#ifndef LIBS_PYCLEO_PYCLEO_ALIASES_HPP_
#define LIBS_PYCLEO_PYCLEO_ALIASES_HPP_

#include <pybind11/pybind11.h>

#include "./optional_terminal_velocity.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_motion.hpp"
#include "cartesiandomain/movement/cartesian_transport_across_domain.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "gridboxes/predcorrmotion.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/gbxindex_observer.hpp"
#include "observers/observers.hpp"
#include "observers/time_observer.hpp"
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
using time = ConstTstepObserver<DoTimeObs<SimpleDataset<FSStore>, FSStore>>;
using gbx = GbxindexObserver<SimpleDataset<FSStore>, FSStore>;
using nullmo = NullSDMMonitor;

using mo = CombinedSDMMonitor<nullmo, nullmo>;
using obs = CombinedObserver<gbx, time, mo>;
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
namespace pycleo_aliases {
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
}  // namespace pycleo_aliases

#endif  // LIBS_PYCLEO_PYCLEO_ALIASES_HPP_
