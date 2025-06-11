/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_observers.cpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 11th June 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for creating python bindings to various different CLEO's Observers instantiations
 */

#include "./py_observers.hpp"

void pyNullObserver(py::module &m) {
  py::class_<pyca::obs_null>(m, "NullObserver")
      .def(py::init())
      .def("next_obs", &pyca::obs_null::next_obs, py::arg("t_mdl"));
}
