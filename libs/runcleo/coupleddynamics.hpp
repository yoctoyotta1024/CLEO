/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: coupleddynamics.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept defining types which can be used for the
 * dynamics solver that couples to SDM in CLEO
 */

#ifndef LIBS_RUNCLEO_COUPLEDDYNAMICS_HPP_
#define LIBS_RUNCLEO_COUPLEDDYNAMICS_HPP_

#include <concepts>

/**
 * @concept CoupledDynamics
 * Concept representing types that meet requirements for dynamics solver coupled to
 * SDM i.e. "CoupledDynamics".
 *
 * A type satisfies the CoupledDynamics concept if it provides the following functions:
 * - prepare_to_timestep(): Prepares for the next timestep.
 * - get_couplstep(): Retrieves the coupling step (unsigned int).
 * - on_step(t): Performs actions associated with a timestep (unsigned int) 't'.
 * - run_step(start, end): Runs the (unsigned int) timestep from 't=start' to 't=end'.
 *
 * @tparam CD The type to check against the CoupledDynamics concept.
 */
template <typename CD>
concept CoupledDynamics = requires(CD cd, unsigned int t) {
  { cd.prepare_to_timestep() } -> std::same_as<void>;
  { cd.get_couplstep() } -> std::convertible_to<unsigned int>;
  { cd.on_step(t) } -> std::same_as<bool>;
  { cd.run_step(t, t) } -> std::same_as<void>;
};

#endif  // LIBS_RUNCLEO_COUPLEDDYNAMICS_HPP_
