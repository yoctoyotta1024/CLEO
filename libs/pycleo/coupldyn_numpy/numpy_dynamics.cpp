/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: numpy_dynamics.cpp
 * Project: coupldyn_numpy
 * Created Date: Wednesday 11th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for struct obeying coupleddynamics concept for dynamics solver in CLEO for
 * coupling between numpy arrays and SDM
 */

#include "./numpy_dynamics.hpp"

void NumpyDynamics::print_dynamics(const unsigned int t_mdl) const {
  for (unsigned int ii = 0; ii < press.size(); ++ii) {
    std::cout << "t: [p, T, qv, qc] = " << t_mdl << ": " << press.data()[ii] << ", "
              << temp.data()[ii] << ", " << qvap.data()[ii] << ", " << qcond.data()[ii] << "\n";
  }
}

void pyNumpyDynamics(py::module &m) {
  py::class_<NumpyDynamics>(m, "NumpyDynamics")
      .def(py::init<const unsigned int, py::array_t<double>, py::array_t<double>,
                    py::array_t<double>, py::array_t<double>, py::array_t<double>,
                    py::array_t<double>, py::array_t<double>>())
      .def_readwrite("press", &NumpyDynamics::press)
      .def_readwrite("temp", &NumpyDynamics::temp)
      .def_readwrite("qvap", &NumpyDynamics::qvap)
      .def_readwrite("qcond", &NumpyDynamics::qcond)
      .def_readwrite("wvel", &NumpyDynamics::wvel)
      .def_readwrite("uvel", &NumpyDynamics::uvel)
      .def_readwrite("vvel", &NumpyDynamics::vvel)
      .def("run_step", &NumpyDynamics::run_step, py::arg("t_mdl"), py::arg("t_next"));
}
