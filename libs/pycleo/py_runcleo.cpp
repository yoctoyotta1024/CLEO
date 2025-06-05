/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_runcleo.cpp
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
 * Functionality for creating python bindings to parts of the runcleo library e.g.
 * to various different CLEO's SDMMethods instantiations
 */

#include "./py_runcleo.hpp"

void pySDMMethods_cartesian_null(py::module &m) {
  py::class_<pyca::CartesianNullSDMMethods>(m, "CartesianNullSDMMethods")
      .def(py::init<const unsigned int, pyca::map_cart, pyca::micro_null, pyca::move_cart_null,
                    pyca::obs_null>())
      .def("get_couplstep", &pyca::CartesianNullSDMMethods::get_couplstep);
}
