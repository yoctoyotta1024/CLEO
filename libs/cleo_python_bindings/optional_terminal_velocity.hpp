/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: optional_terminal_velocity.hpp
 * Project: cleo_python_bindings
 * Created Date: Friday 11th July 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_OPTIONAL_TERMINAL_VELOCITY_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_OPTIONAL_TERMINAL_VELOCITY_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"

/**
 * @brief Different options for terminal velocity combined into one struct.
 */
class OptionalTerminalVelocity {
 private:
  bool enable_terminal_velocity;
  RogersGKTerminalVelocity rogersgk;

 public:
  explicit OptionalTerminalVelocity(const bool enable_terminal_velocity)
      : enable_terminal_velocity(enable_terminal_velocity), rogersgk(RogersGKTerminalVelocity{}) {}

  /**
   * @brief Returns the terminal velocity of a droplet from various different options.
   *
   * @param drop The superdroplet.
   * @return The (dimensionless) terminal velocity.
   */
  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop) const {
    if (enable_terminal_velocity) {
      return rogersgk(drop);
    }
    return 0.0;
  }
};

#endif  // LIBS_CLEO_PYTHON_BINDINGS_OPTIONAL_TERMINAL_VELOCITY_HPP_
