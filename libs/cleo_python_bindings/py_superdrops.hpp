/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_superdrops.hpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Python bindings to various parts of CLEO's superdrops library
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_PY_SUPERDROPS_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_PY_SUPERDROPS_HPP_

#include <pybind11/pybind11.h>

#include <iostream>
#include <limits>

#include "../cleoconstants.hpp"
#include "./cleo_python_bindings_aliases.hpp"
#include "configuration/config.hpp"
#include "initialise/timesteps.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/collisions.hpp"
#include "superdrops/collisions/longhydroprob.hpp"
#include "superdrops/condensation.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"

namespace py = pybind11;
namespace pyca = cleo_python_bindings_aliases;

void pyNullMicrophysicalProcess(py::module &m);
void pyAllMicrophysicalProcess(py::module &m);
void pycreate_microphysical_process(py::module &m);

void pyNullMotion(py::module &m);

#endif  // LIBS_CLEO_PYTHON_BINDINGS_PY_SUPERDROPS_HPP_
