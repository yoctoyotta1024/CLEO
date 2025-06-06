/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_superdrops.cpp
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
 * Functionality for creating python bindings to various parts of CLEO's superdrops library
 */

#include "./py_superdrops.hpp"

void pyMicrophysicalProcess_null(py::module &m) {
  py::class_<pyca::micro_null>(m, "NullMicrophysicalProcess")
      .def(py::init());
}

void pyMotion_null(py::module &m) {
  py::class_<pyca::mo_null>(m, "NullMotion")
      .def(py::init());
}
