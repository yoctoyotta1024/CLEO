/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_initialise.hpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 10th June 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Python bindings to various different CLEO initialisation/configuration functions and structures
 */

#ifndef LIBS_PYCLEO_PY_INITIALISE_HPP_
#define LIBS_PYCLEO_PY_INITIALISE_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>

#include "initialise/config.hpp"
#include "initialise/init_supers_from_binary.hpp"
#include "initialise/initgbxsnull.hpp"
#include "initialise/optional_config_params.hpp"
#include "initialise/timesteps.hpp"

namespace py = pybind11;

void pyConfig(py::module &m);

void pyTimesteps(py::module &m);
void pycreate_timesteps(py::module &m);

void pyInitSupersFromBinary(py::module &m);
void pyInitSupersFromBinaryParams(py::module &m);

void pyInitGbxsNull(py::module &m);

#endif  // LIBS_PYCLEO_PY_INITIALISE_HPP_
