/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cvodecomms.hpp
 * Project: coupldyn_cvode
 * Created Date: Friday 27th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * send and receive dynamics functions for SDM when coupled to the CVODE ODE solver
 */

#ifndef LIBS_COUPLDYN_CVODE_CVODECOMMS_HPP_
#define LIBS_COUPLDYN_CVODE_CVODECOMMS_HPP_

#include <array>
#include <vector>

#include "../kokkosaliases.hpp"
#include "coupldyn_cvode/cvodedynamics.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"

struct CvodeComms {
 private:
  /* get change in state since
  previous time step to current one */
  std::array<double, 4> state_change(CvodeDynamics &cvode, const viewh_constgbx h_gbxs,
                                     const size_t ii) const;

  /* change is_delta_y = false to is_delta_y = true
  if delta contains non-zero elements */
  bool is_state_change(const std::array<double, 4> &delta, bool is_delta_y) const;

 public:
  /*
  update Gridboxes' states using information received from CVODE dynanmics
  solver for  press, temp, qvap and qcond.
  Note: ii indexing for cvode isn't compatible with MPI domain decompositon.
  */
  template <typename GbxMaps, typename CD = CvodeDynamics>
  void receive_dynamics(const GbxMaps &gbxmaps, const CvodeDynamics &cvode,
                        const viewh_gbx h_gbxs) const {
    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii) {
      State &state(h_gbxs(ii).state);
      const auto cvodestate(cvode.get_current_state(ii));  // ii'th states' [p, t, qv, qc]

      state.press = cvodestate.at(0);
      state.temp = cvodestate.at(1);
      state.qvap = cvodestate.at(2);
      state.qcond = cvodestate.at(3);
    }
  }

  /*
  send information from Gridboxes' states to CVODE dynanmics solver for  temp, qvap
  and qcond (excludes press)
  Note: ii indexing for cvode isn't compatible with MPI domain decompositon.
  */
  template <typename GbxMaps, typename CD = CvodeDynamics>
  void send_dynamics(const GbxMaps &gbxmaps, const viewh_constgbx h_gbxs,
                     CvodeDynamics &cvode) const {
    std::vector<double> delta_y;
    bool is_delta_y(false);

    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii) {
      const auto delta(state_change(cvode, h_gbxs, ii));  // ii'th [press, temp, qvap, qcond] change

      delta_y.push_back(delta.at(0));
      delta_y.push_back(delta.at(1));
      delta_y.push_back(delta.at(2));
      delta_y.push_back(delta.at(3));

      is_delta_y = is_state_change(delta, is_delta_y);
    }

    if (is_delta_y) {
      cvode.reinitialise(cvode.get_time(), delta_y);
    }
  }
};

#endif  // LIBS_COUPLDYN_CVODE_CVODECOMMS_HPP_
