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
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Entry point for CLEO's python bindngs for the coupldyn_numpy library
 * sub-module of cleo_python_bindings
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_PYCOUPLDYN_NUMPY_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_PYCOUPLDYN_NUMPY_HPP_

#include <pybind11/pybind11.h>

#include <iostream>

#include "./numpy_comms.hpp"
#include "./numpy_dynamics.hpp"

namespace py = pybind11;

int test_coupldyn_numpy(const int i, const int j);

void include_coupldyn_numpy_submodule(py::module &m);

#endif  // LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_PYCOUPLDYN_NUMPY_HPP_
