/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_configuration.cpp
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
 * Functionality for creating python bindings to various different CLEO
 * configuration functions and structures
 */

#include "./py_configuration.hpp"

void pyConfig(py::module &m) {
  py::class_<Config>(m, "Config")
      .def(py::init<const std::filesystem::path>(), py::arg("config_filename"))
      .def("get_ngbxs", &Config::get_ngbxs)
      .def("get_nspacedims", &Config::get_nspacedims)
      .def("get_grid_filename", &Config::get_grid_filename)
      .def("get_initsupersfrombinary", &Config::get_initsupersfrombinary);
}

void pyInitSupersFromBinaryParams(py::module &m) {
  py::class_<OptionalConfigParams::InitSupersFromBinaryParams>(m, "InitSupersFromBinaryParams");
}
