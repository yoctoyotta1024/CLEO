/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_initialise.hpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Python bindings to various different CLEO initialisation functions and structures
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_PY_INITIALISE_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_PY_INITIALISE_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>

#include "configuration/config.hpp"
#include "initialise/init_supers_from_binary.hpp"
#include "initialise/initgbxsnull.hpp"
#include "initialise/timesteps.hpp"

namespace py = pybind11;

void pyTimesteps(py::module &m);
void pycreate_timesteps(py::module &m);
void pyrealtime2step(py::module &m);

void pyInitSupersFromBinary(py::module &m);

void pyInitGbxsNull(py::module &m);

#endif  // LIBS_CLEO_PYTHON_BINDINGS_PY_INITIALISE_HPP_
