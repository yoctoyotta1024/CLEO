/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_cartesiandomain.hpp
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
 * Python bindings to various parts of CLEO's cartesiandomain library
 */

#ifndef LIBS_PYCLEO_PY_CARTESIANDOMAIN_HPP_
#define LIBS_PYCLEO_PY_CARTESIANDOMAIN_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>

#include "./pycleo_aliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"

namespace py = pybind11;
namespace pyca = pycleo_aliases;

void pyCartesianMaps(py::module &m);

void pycreate_cartesian_maps(py::module &m);

#endif  // LIBS_PYCLEO_PY_CARTESIANDOMAIN_HPP_
