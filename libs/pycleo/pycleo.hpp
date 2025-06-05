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
 * Last Modified: Thursday 5th June 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Entry point for creating CLEO's python bindngs
 */

#ifndef LIBS_PYCLEO_PYCLEO_HPP_
#define LIBS_PYCLEO_PYCLEO_HPP_

#include <pybind11/pybind11.h>

#include <Kokkos_Core.hpp>
#include <iostream>

namespace py = pybind11;

#include "runcleo/sdmmethods.hpp"

int test_python_bindings(const int i, const int j);

PYBIND11_MODULE(pycleo, m) {
  m.doc() = "Python bindings for selected parts of CLEO's libraries";

  m.def("test_python_bindings", &test_python_bindings, "test function for CLEO example",
        py::arg("i"), py::arg("j"));
}

#endif  // LIBS_PYCLEO_PYCLEO_HPP_
