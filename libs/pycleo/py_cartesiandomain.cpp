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
 * Last Modified: Friday 6th June 2025
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
  py::class_<pyca::map_cart>(m, "CartesianMaps").def(py::init());
}

void pycreate_cartesian_maps(py::module &m) {
  m.def("create_cartesian_maps", &create_cartesian_maps, "returns CartesianMaps instance",
        py::arg("ngbxs"), py::arg("nspacedims"), py::arg("grid_filename"));
}

void pyCartesianTransportAcrossDomain(py::module &m) {
  py::class_<pyca::trans_cart>(m, "CartesianTransportAcrossDomain").def(py::init());
}
