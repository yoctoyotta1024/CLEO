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
 * Last Modified: Tuesday 1st July 2025
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

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_motion.hpp"
#include "cartesiandomain/movement/cartesian_transport_across_domain.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "gridboxes/predcorrmotion.hpp"
#include "observers/observers.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/condensation.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/terminalvelocity.hpp"

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
using micro_all =
    CombinedMicrophysicalProcess<NullMicrophysicalProcess, ConstTstepMicrophysics<DoCondensation>>;

using mo_null = NullMotion;
using mo_cart_predcorr =
    PredCorrMotion<CartesianMaps, RogersGKTerminalVelocity, CartesianCheckBounds>;
using bcs_null = NullBoundaryConditions;
using trans_cart = CartesianTransportAcrossDomain;
using move_cart_null = MoveSupersInDomain<map_cart, mo_null, trans_cart, bcs_null>;
using move_cart = MoveSupersInDomain<map_cart, mo_cart_predcorr, trans_cart, bcs_null>;

using obs_null = NullObserver;

using sdm_cart_null = SDMMethods<map_cart, micro_null, mo_null, trans_cart, bcs_null, obs_null>;
using sdm_cart_all =
    SDMMethods<map_cart, micro_all, mo_cart_predcorr, trans_cart, bcs_null, obs_null>;
}  // namespace pycleo_aliases

#endif  // LIBS_PYCLEO_PYCLEO_ALIASES_HPP_
