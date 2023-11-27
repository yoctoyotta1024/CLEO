/*
 * ----- CLEO -----
 * File: fromfilecomms.hpp
 * Project: coupldyn_fromfile
 * Created Date: Tuesday 31st October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 7th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * send and receive dynamics functions
 * for SDM when coupled to the fromfile
 * dynamics solver
 */

#ifndef FROMFILECOMMS_HPP
#define FROMFILECOMMS_HPP

#include "../kokkosaliases.hpp"
#include "./fromfile_cartesian_dynamics.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"

struct FromFileComms
/* 1-way coupling from coupldyn to CLEO's gridboxes where
coupldyn sends information to h_gbxs but doesn't
receive any back. Struct obeys coupling comms concept */
{
private:
  void update_gridbox_state(const FromFileDynamics &ffdyn,
                            const size_t ii,
                            Gridbox &gbx) const;
  /* updates the state of a gridbox using information
  received from FromFileDynamics solver for 1-way
  coupling to CLEO SDM */ 

public:
  template <typename CD = FromFileDynamics>
  void send_dynamics(const viewh_constgbx h_gbxs,
                     FromFileDynamics &ffdyn) const {}
  /* send information from Gridboxes' states
  to coupldyn is null for FromFileDynamics*/

  template <typename CD = FromFileDynamics>
  void receive_dynamics(const FromFileDynamics &ffdyn,
                        const viewh_gbx h_gbxs) const;
  /* update Gridboxes' states using information
  received from FromFileDynamics solver for
  1-way coupling to CLEO SDM */


};

#endif // FROMFILECOMMS_HPP