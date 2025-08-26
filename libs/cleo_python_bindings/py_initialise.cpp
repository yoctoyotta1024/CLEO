/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_initialise.cpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for creating python bindings to various different CLEO
 * initialisation functions and structures
 */

#include "./py_initialise.hpp"

void pyTimesteps(py::module &m) {
  py::class_<Timesteps>(m, "Timesteps")
      .def(py::init<const RequiredConfigParams::TimestepsParams &>())
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

void pyrealtime2step(py::module &m) {
  m.def(
      "realtime2step", [](const double TSTEP) { return realtime2step(TSTEP); },
      "converts time [s] into SDM model timestep", py::arg("TSTEP"));
}

void pyInitSupersFromBinary(py::module &m) {
  py::class_<InitSupersFromBinary>(m, "InitSupersFromBinary")
      .def(py::init<const OptionalConfigParams::InitSupersFromBinaryParams &,
                    const CartesianMaps &>());
}

void pyInitGbxsNull(py::module &m) {
  py::class_<InitGbxsNull>(m, "InitGbxsNull").def(py::init<const size_t>());
}
