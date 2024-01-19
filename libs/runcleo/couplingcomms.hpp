/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: couplingcomms.hpp
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
 * File Description:
 * concept defining types which can be used for the
 * coupling between SDMMethods and a dynamics solver
 * in RunCLEO
 */


#ifndef COUPLINGCOMMS_HPP
#define COUPLINGCOMMS_HPP

#include <concepts>

#include "../kokkosaliases.hpp"
#include "./coupleddynamics.hpp"

template <typename Comms, typename CD>
concept CouplingComms = requires(Comms s,
                            CD &coupldyn,
                            viewh_gbx h_gbxs) {
  {
    s.template send_dynamics<CD>(h_gbxs, coupldyn)
  } -> std::same_as<void>;
  {
    s.template receive_dynamics<CD>(coupldyn, h_gbxs)
  } -> std::same_as<void>;
};

struct NullComms
/* null coupling doesnt send or receive information between
coupldyn and h_gbxs but still obeys coupling comms concept */
{
  template <CoupledDynamics CD>
  void receive_dynamics(const CD &coupldyn,
                        const viewh_gbx h_gbxs) const {}
  /* update Gridboxes' states using information
  received from coupldyn */

  template <CoupledDynamics CD>
  void send_dynamics(const viewh_constgbx h_gbxs,
                     CD &coupldyn) const {}
  /* send information from Gridboxes' states to coupldyn */
};

#endif // COUPLINGCOMMS_HPP
