/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_initialise.cpp
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
 * Functionality for creating python bindings to various different CLEO
 * initialisation/configuration functions and structures
 */

#include "./py_initialise.hpp"

void pyConfig(py::module &m) {
  py::class_<Config>(m, "Config")
      .def(py::init<const std::filesystem::path>(), py::arg("config_filename"));
}

void pyTimesteps(py::module &m) {
  py::class_<Timesteps>(m, "Timesteps")
      .def(py::init<const RequiredConfigParams::TimestepsParams>())
      .def("get_condstep", &Timesteps::get_condstep)
      .def("get_collstep", &Timesteps::get_collstep)
      .def("get_motionstep", &Timesteps::get_motionstep)
      .def("get_couplstep", &Timesteps::get_couplstep)
      .def("get_obsstep", &Timesteps::get_obsstep)
      .def("get_t_end", &Timesteps::get_t_end);
}

void pycreate_timesteps(py::module &m) {
  m.def(
      "pycreate_timesteps", [](const Config &config) { return Timesteps(config.get_timesteps()); },
      "returns Timesteps instance", py::arg("config"));
}
