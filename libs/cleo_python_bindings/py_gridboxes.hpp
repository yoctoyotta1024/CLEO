/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_gridboxes.hpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Python bindings to various parts of CLEO's gridboxes library
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_PY_GRIDBOXES_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_PY_GRIDBOXES_HPP_

#include <pybind11/pybind11.h>

#include "../kokkosaliases.hpp"
#include "./cleo_python_bindings_aliases.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "gridboxes/supersindomain.hpp"

namespace py = pybind11;
namespace pyca = cleo_python_bindings_aliases;

void pyNullBoundaryConditions(py::module &m);

void pyCartesianNullMoveSupersInDomain(py::module &m);

void pySupersInDomain(py::module &m);

void pyGridboxesDualView(py::module &m);

#endif  // LIBS_CLEO_PYTHON_BINDINGS_PY_GRIDBOXES_HPP_
