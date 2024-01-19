/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: coupleddynamics.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept defining types which can be used for the
 * dynamics solver that coupled to SDM in CLEO
 */

#ifndef COUPLEDDYNAMICS_HPP
#define COUPLEDDYNAMICS_HPP

#include <concepts>

template <typename CD>
concept CoupledDynamics = requires(CD cd, unsigned int t)
/* concept Coupled Dynamics is all types that meet requirements
(constraints) of preparation and timestepping functions
and have constant unsigned int type (for interval) */
{
  {
    cd.prepare_to_timestep()
  } -> std::same_as<void>;
  {
    cd.get_couplstep()
  } -> std::convertible_to<unsigned int>;
  {
    cd.on_step(t)
  } -> std::same_as<bool>;
  {
    cd.run_step(t, t)
  } -> std::same_as<void>;
};

#endif // COUPLEDDYNAMICS_HPP
