/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: pycoupldyn_numpy.cpp
 * Project: coupldyn_numpy
 * Created Date: Wednesday 11th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for creating CLEO's python bindngs for the coupldyn_numpy
 * library sub-module of cleo_python_bindings
 */

#include "./pycoupldyn_numpy.hpp"

int test_coupldyn_numpy(const int i, const int j) {
  std::cout << "Hello World\n";
  return i * j;
}

void include_coupldyn_numpy_submodule(py::module &m) {
  auto m_sub = m.def_submodule(
      "coupldyn_numpy", "Python bindings for selected parts of CLEO's coupldyn_numpy library");

  m_sub.def("test_coupldyn_numpy", &test_coupldyn_numpy,
            "test function for coupldyn_numpy sub-module", py::arg("i"), py::arg("j"));

  /* dynamics */
  pyNumpyDynamics(m_sub);

  /* coupling */
  pyNumpyComms(m_sub);
}
