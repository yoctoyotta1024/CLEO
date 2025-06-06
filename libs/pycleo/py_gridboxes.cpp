/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_gridboxes.cpp
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
 * Functionality for creating python bindings to various parts of CLEO's gridboxes library
 */

#include "./py_gridboxes.hpp"

void pyBoundaryConditions_null(py::module &m) {
  py::class_<pyca::bcs_null>(m, "NullBoundaryConditions")
      .def(py::init());
}

void pyMoveSupersInDomain_cartesian_null(py::module &m) {
  py::class_<pyca::move_cart_null>(m, "CartesianNullMoveSupersInDomain")
      .def(py::init<pyca::mo_null, pyca::trans_cart, pyca::bcs_null>());
}
