/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_gridboxes.hpp
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
 * Python bindings to various parts of CLEO's gridboxes library
 */

#ifndef LIBS_PYCLEO_PY_GRIDBOXES_HPP_
#define LIBS_PYCLEO_PY_GRIDBOXES_HPP_

#include <pybind11/pybind11.h>

#include "./pycleo_aliases.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/movesupersindomain.hpp"

namespace py = pybind11;
namespace pyca = pycleo_aliases;

void pyBoundaryConditions_null(py::module &m);

void pyMoveSupersInDomain_cartesian_null(py::module &m);

#endif  // LIBS_PYCLEO_PY_GRIDBOXES_HPP_
