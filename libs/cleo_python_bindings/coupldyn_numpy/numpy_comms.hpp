/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: numpy_comms.hpp
 * Project: coupldyn_numpy
 * Created Date: Friday 17th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct obeying coupling concept for dynamics solver in CLEO for
 * coupling between NumpyDynamics and SDM
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_NUMPY_COMMS_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_NUMPY_COMMS_HPP_

#include <Kokkos_Core.hpp>

#include "../../kokkosaliases.hpp"
#include "./numpy_dynamics.hpp"
#include "cartesiandomain/cartesianmaps.hpp"

void pyNumpyComms(py::module &m);

/* coupling of NumpyDynamics to CLEO's gridboxes. Struct obeys coupling comms concept */
struct NumpyComms {
  /* receive information from NumpyDynamics solver to CLEO SDM */
  template <typename GbxMaps, typename CD = NumpyComms>
  KOKKOS_FUNCTION void receive_dynamics(const GbxMaps &gbxmaps, const NumpyDynamics &numpydyn,
                                        const viewh_gbx h_gbxs) const;

  /* send information from Gridboxes' states to NumpyDynamics */
  template <typename GbxMaps, typename CD = NumpyComms>
  KOKKOS_FUNCTION void send_dynamics(const GbxMaps &gbxmaps, const viewh_constgbx h_gbxs,
                                     NumpyDynamics &numpydyn) const;
};

#endif  // LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_NUMPY_COMMS_HPP_
