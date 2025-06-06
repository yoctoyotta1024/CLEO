/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: pycleo.hpp
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
 * Entry point for creating CLEO's python bindngs module
 */

#ifndef LIBS_PYCLEO_PYCLEO_HPP_
#define LIBS_PYCLEO_PYCLEO_HPP_

#include <pybind11/pybind11.h>

#include <iostream>

#include "./py_observers.hpp"
#include "./py_runcleo.hpp"
#include "./py_cartesiandomain.hpp"
// #include "./py_gridboxes.hpp" #WIP
// #include "./py_superdrops.hpp" #WIP

namespace py = pybind11;

int test_python_bindings(const int i, const int j);

PYBIND11_MODULE(pycleo, m) {
  m.doc() = "Python bindings for selected parts of CLEO's libraries";

  m.def("test_python_bindings", &test_python_bindings, "test function for CLEO example",
        py::arg("i"), py::arg("j"));

  /* maps */
  pyCartesianMaps(m);
  pycreate_cartesian_maps(m);

  /* microphyiscs */
  // pyMicrophysicalProcess_null(m); #WIP

  /* motion */
  // pyMotion_null(m); #WIP

  /* boundary conditions */
  // pyBoundaryConditions_null(m); #WIP

  /* transport */
  // pyTransportAcrossDomain_cartesian(m); #WIP

  /* movement */
  // pyMoveSupersInDomain_cartesian_null(m); #WIP

  /* observers */
  pyObservers_null(m);

  /* sdmmethods */
  pySDMMethods_cartesian_null(m);
}

#endif  // LIBS_PYCLEO_PYCLEO_HPP_
