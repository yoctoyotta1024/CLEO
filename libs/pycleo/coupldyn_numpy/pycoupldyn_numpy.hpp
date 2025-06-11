/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: pycoupldyn_numpy.hpp
 * Project: coupldyn_numpy
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
 * Entry point for CLEO's python bindngs module for the coupldyn_numpy library
 */

#ifndef LIBS_PYCLEO_COUPLDYN_NUMPY_PYCOUPLDYN_NUMPY_HPP_
#define LIBS_PYCLEO_COUPLDYN_NUMPY_PYCOUPLDYN_NUMPY_HPP_

#include <mpi.h>
#include <pybind11/pybind11.h>

#include <iostream>

#include "./numpy_comms.hpp"
#include "./numpy_dynamics.hpp"

namespace py = pybind11;

int test_python_bindings(const int i, const int j);

PYBIND11_MODULE(coupldyn_numpy, m) {
  m.doc() = "Python bindings for selected parts of CLEO's coupldyn_numpy library";

  m.def("test_python_bindings", &test_python_bindings, "test function for coupldyn_numpy module",
        py::arg("i"), py::arg("j"));

  /* dynamics */
  pyNumpyDynamics(m);

  /* coupling */
  pyNumpyComms(m);
}

#endif  // LIBS_PYCLEO_COUPLDYN_NUMPY_PYCOUPLDYN_NUMPY_HPP_
