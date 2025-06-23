/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_configuration.hpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 12th June 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Python bindings to various different CLEO configuration functions and structures
 */

#ifndef LIBS_PYCLEO_PY_CONFIGURATION_HPP_
#define LIBS_PYCLEO_PY_CONFIGURATION_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>

#include "configuration/config.hpp"
#include "configuration/optional_config_params.hpp"

namespace py = pybind11;

void pyConfig(py::module &m);

void pyInitSupersFromBinaryParams(py::module &m);

#endif  // LIBS_PYCLEO_PY_CONFIGURATION_HPP_
