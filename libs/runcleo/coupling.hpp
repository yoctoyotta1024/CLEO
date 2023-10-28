/*
 * ----- CLEO -----
 * File: coupling.hpp
 * Project: runcleo
 * Created Date: Sunday 29th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 29th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * concept defining types which can be used for the
 * coupling between SDMMethods and a dynamics solver 
 * in RunCLEO
 */


#ifndef COUPLING_HPP 
#define COUPLING_HPP 

#include <concepts>

#include "../kokkosaliases.hpp"

// TODO finish this concept with templated type CD

// template <typename Comms, CoupledDynamics CD>
// concept Coupling = requires(Comms comms, CD &coupldyn,
//                             viewh_gbx h_gbxs)
// /* concept Coupling is all types that meet requirements
// (constraints) of communicaiton functions */
// {
//   {
//     comms.send_dynamics(h_gbxs, coupldyn)
//   } -> std::same_as<void>;
//   {
//     comms.receive_dynamics(coupldyn, h_gbxs)
//   } -> std::same_as<void>;
// };

template <CoupledDynamics CD>
struct NullComms
{
  NullComms(const CD coupldyn) {}

  void receive_dynamics(const CD &coupldyn,
                        const viewh_gbx h_gbxs) const {}
  /* update Gridboxes' states using information
  received from coupldyn */

  void send_dynamics(const viewh_constgbx h_gbxs,
                     CD &coupldyn) const {}
  /* send information from Gridboxes' states to coupldyn */
};

#endif // COUPLING_HPP  