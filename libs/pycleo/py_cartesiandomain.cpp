/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_cartesiandomain.cpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 1st July 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for creating python bindings to various parts of CLEO's cartesiandomain library
 */

#include "./py_cartesiandomain.hpp"

void pyCartesianMaps(py::module &m) {
  py::class_<pyca::map_cart>(m, "CartesianMaps")
      .def(py::init())
      .def("get_local_ngridboxes_hostcopy", &pyca::map_cart::get_local_ngridboxes_hostcopy);
}

void pycreate_cartesian_maps(py::module &m) {
  m.def("create_cartesian_maps", &create_cartesian_maps, "returns CartesianMaps instance",
        py::arg("ngbxs"), py::arg("nspacedims"), py::arg("grid_filename"));
}

void pyCartesianTransportAcrossDomain(py::module &m) {
  py::class_<pyca::trans_cart>(m, "CartesianTransportAcrossDomain").def(py::init());
}

void pyCartesianPredCorrMotion(py::module &m) {
  py::class_<pyca::mo_cart_predcorr>(m, "CartesianPredCorrMoveSupersInDomain")
      .def(py::init<unsigned int, std::function<double(unsigned int)>, RogersGKTerminalVelocity,
                    CartesianCheckBounds>());
}

void pycreate_cartesian_predcorr_motion(py::module &m) {
  m.def(
      "create_cartesian_predcorr_motion",
      [](const unsigned int motionstep) {
        return CartesianMotion(motionstep, &step2dimlesstime, RogersGKTerminalVelocity{});
      },
      "returns CartesianPredCorrMotion instance", py::arg("motionstep"));
}

void pyCartesianMoveSupersInDomain(py::module &m) {
  py::class_<pyca::move_cart>(m, "CartesianMoveSupersInDomain")
      .def(py::init<pyca::mo_cart_predcorr, pyca::trans_cart, pyca::bcs_null>());
}
