/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: couplingcomms.hpp
 * Project: runcleo
 * Created Date: Sunday 29th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept defining types which can be used for the
 * coupling between SDMMethods and a dynamics solver
 * in RunCLEO, also null instance for coupling called
 * "NullComms"
 */

#ifndef LIBS_RUNCLEO_COUPLINGCOMMS_HPP_
#define LIBS_RUNCLEO_COUPLINGCOMMS_HPP_

#include <concepts>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "runcleo/coupleddynamics.hpp"

/**
 * @concept CouplingComms
 * Concept representing types that handle communication between SDM and coupled dynamics solver.
 *
 * A type satisfies the CouplingComms concept if it provides the following functions:
 * - `send_dynamics(h_gbxs, coupldyn)`: Sends dynamics information from SDM view of Gridboxes
 * `h_gbxs` to coupled dynamics solver `coupldyn`.
 * - `receive_dynamics(coupldyn, h_gbxs)`: Receives dynamics information from coupled dynamics
 * solver `coupldyn` into SDM view of Gridboxes `h_gbxs`.
 *
 * @tparam Comms The type for communication to check against the CouplingComms concept.
 * @tparam GbxMaps The type for gridbox maps to check against the GridboxMaps concept.
 * @tparam CD The type for the dyanmics solver to check against the CoupledDynamics concept.
 */
template <typename Comms, typename GbxMaps, typename CD>
concept CouplingComms = requires(Comms s, GbxMaps &gbxmaps, CD &coupldyn, viewh_gbx h_gbxs) {
  { s.template send_dynamics<GbxMaps, CD>(gbxmaps, h_gbxs, coupldyn) } -> std::same_as<void>;
  { s.template receive_dynamics<GbxMaps, CD>(gbxmaps, coupldyn, h_gbxs) } -> std::same_as<void>;
};

/**
 * @struct NullComms
 * @brief Represents a null communication handler that doesn't send or receive information.
 *
 * The NullComms struct implements the CouplingComms concept but doesn't perform any communication
 * between SDM Gridboxes and CoupledDynamics solver.
 */
struct NullComms {
  /**
   * @brief Receives dynamics information.
   *
   * This function does nothing as it represents a null communication handler.
   *
   * @tparam CD The coupled dynamics solver type.
   * @param gbxmaps The Gridbox Maps.
   * @param coupldyn The coupled dynamics solver object.
   * @param h_gbxs The view of Gridboxes.
   */
  template <GridboxMaps GbxMaps, CoupledDynamics CD>
  void receive_dynamics(const GbxMaps &gbxmaps, const CD &coupldyn, const viewh_gbx h_gbxs) const {}

  /**
   * @brief Sends dynamics information.
   *
   * This function does nothing as it represents a null communication handler.
   *
   * @tparam CD The coupled dynamics solver type.
   * @param gbxmaps The Gridbox Maps.
   * @param h_gbxs The view of Gridboxes.
   * @param coupldyn The coupled dynamics solver object.
   */
  template <GridboxMaps GbxMaps, CoupledDynamics CD>
  void send_dynamics(const GbxMaps &gbxmaps, const viewh_constgbx h_gbxs, CD &coupldyn) const {}
};

#endif  // LIBS_RUNCLEO_COUPLINGCOMMS_HPP_
