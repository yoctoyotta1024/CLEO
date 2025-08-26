/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_zarr.hpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Python bindings to parts of CLEO's zarr library
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_PY_ZARR_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_PY_ZARR_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>

#include "./cleo_python_bindings_aliases.hpp"
#include "configuration/config.hpp"
#include "zarr/fsstore.hpp"
#include "zarr/simple_dataset.hpp"

namespace py = pybind11;
namespace pyca = cleo_python_bindings_aliases;

inline void pyFSStore(py::module &m) {
  py::class_<FSStore>(m, "FSStore").def(py::init<std::filesystem::path>());
}

inline void pySimpleDataset(py::module &m) {
  py::class_<SimpleDataset<FSStore>>(m, "SimpleDataset").def(py::init<FSStore &>());
}

#endif  // LIBS_CLEO_PYTHON_BINDINGS_PY_ZARR_HPP_
