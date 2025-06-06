/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_observers.hpp
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
 * Python bindings to various different CLEO's Observers instantiations
 */

#ifndef LIBS_PYCLEO_PY_OBSERVERS_HPP_
#define LIBS_PYCLEO_PY_OBSERVERS_HPP_

#include <pybind11/pybind11.h>

#include "./pycleo_aliases.hpp"
#include "observers/observers.hpp"

namespace py = pybind11;
namespace pyca = pycleo_aliases;

void pyNullObserver(py::module &m);

#endif  // LIBS_PYCLEO_PY_OBSERVERS_HPP_
