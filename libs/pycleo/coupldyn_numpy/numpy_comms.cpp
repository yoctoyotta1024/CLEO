/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: numpy_comms.cpp
 * Project: coupldyn_numpy
 * Created Date: Wednesday 11th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for struct obeying coupling concept for dynamics solver in CLEO for
 * coupling between NumpyDynamics and SDM
 */

#include "./numpy_comms.hpp"

/* receive information from NumpyDynamics solver to CLEO SDM */
template <typename GbxMaps, typename CD>
KOKKOS_FUNCTION void NumpyComms::receive_dynamics(const GbxMaps &gbxmaps,
                                                  const NumpyDynamics &numpydyn,
                                                  const viewh_gbx h_gbxs) const {
  const size_t ngbxs(h_gbxs.extent(0));

  Kokkos::parallel_for("receive_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs),
                       [=](const size_t ii) {
                         // for (size_t ii = 0; ii < ngbxs; ++ii) {
                         const auto idx = gbxmaps.local_to_global_gridbox_index(ii);
                         State &state(h_gbxs(ii).state);

                         state.press = numpydyn.get_press(idx);
                         state.temp = numpydyn.get_temp(idx);
                         state.qvap = numpydyn.get_qvap(idx);
                         state.qcond = numpydyn.get_qcond(idx);

                         state.wvel = numpydyn.get_wvel(ii);
                         state.uvel = numpydyn.get_uvel(ii);
                         state.vvel = numpydyn.get_vvel(ii);
                       });
}

/* send information from Gridboxes' states to NumpyDynamics */
template <typename GbxMaps, typename CD>
KOKKOS_FUNCTION void NumpyComms::send_dynamics(const GbxMaps &gbxmaps, const viewh_constgbx h_gbxs,
                                               NumpyDynamics &numpydyn) const {
  const size_t ngbxs(h_gbxs.extent(0));

  Kokkos::parallel_for("send_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs),
                       [=, &numpydyn](const size_t ii) {
                         const auto idx = gbxmaps.local_to_global_gridbox_index(ii);
                         State &state(h_gbxs(ii).state);

                         numpydyn.set_press(idx, state.press);
                         numpydyn.set_temp(idx, state.temp);
                         numpydyn.set_qvap(idx, state.qvap);
                         numpydyn.set_qcond(idx, state.qcond);
                       });
}

template void NumpyComms::receive_dynamics<CartesianMaps, NumpyComms>(const CartesianMaps &,
                                                                      const NumpyDynamics &,
                                                                      const viewh_gbx) const;

template void NumpyComms::send_dynamics<CartesianMaps, NumpyComms>(const CartesianMaps &,
                                                                   const viewh_constgbx,
                                                                   NumpyDynamics &) const;

void pyNumpyComms(py::module &m) {
  py::class_<NumpyComms>(m, "NumpyComms")
      .def(py::init())
      .def(
          "receive_dynamics",
          [](const NumpyComms &self, const CartesianMaps &gbxmaps, const NumpyDynamics &numpydyn,
             const dualview_gbx gbxs) {
            self.receive_dynamics(gbxmaps, numpydyn, gbxs.view_host());
          },
          py::arg("gbxmaps"), py::arg("numpydyn"), py::arg("h_gbxs"))
      .def(
          "send_dynamics",
          [](const NumpyComms &self, const CartesianMaps &gbxmaps, const dualview_gbx gbxs,
             NumpyDynamics &numpydyn) { self.send_dynamics(gbxmaps, gbxs.view_host(), numpydyn); },
          py::arg("gbxmaps"), py::arg("h_gbxs"), py::arg("numpydyn"));
}
