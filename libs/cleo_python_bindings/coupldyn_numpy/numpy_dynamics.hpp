/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: numpy_dynamics.hpp
 * Project: coupldyn_numpy
 * Created Date: Wednesday 11th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct obeying coupleddynamics concept for dynamics solver in CLEO for
 * coupling betaween numpy arrays and SDM
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_NUMPY_DYNAMICS_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_NUMPY_DYNAMICS_HPP_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <iostream>
#include <utility>

namespace py = pybind11;

void pyNumpyDynamics(py::module &m);

struct NumpyDynamics {
 private:
  const unsigned int interval;
  void print_dynamics(const unsigned int t_mdl) const;

 public:
  py::array_t<double> press;
  py::array_t<double> temp;
  py::array_t<double> qvap;
  py::array_t<double> qcond;
  py::array_t<double> wvel;
  py::array_t<double> uvel;
  py::array_t<double> vvel;

  explicit NumpyDynamics(const unsigned int couplstep, py::array_t<double> press,
                         py::array_t<double> temp, py::array_t<double> qvap,
                         py::array_t<double> qcond, py::array_t<double> wvel,
                         py::array_t<double> uvel, py::array_t<double> vvel)
      : interval(couplstep),
        press(press),
        temp(temp),
        qvap(qvap),
        qcond(qcond),
        wvel(wvel),
        uvel(uvel),
        vvel(vvel) {}

  void prepare_to_timestep() const {}

  unsigned int get_couplstep() const { return interval; }

  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  void run_step(const unsigned int t_mdl, const unsigned int t_next) const {
    if (on_step(t_mdl)) {
      // print_dynamics(t_mdl);  // useful for debugging
    }
  }

  double get_press(const size_t ii) const { return press.data()[ii]; }

  double get_temp(const size_t ii) const { return temp.data()[ii]; }

  double get_qvap(const size_t ii) const { return qvap.data()[ii]; }

  double get_qcond(const size_t ii) const { return qcond.data()[ii]; }

  std::pair<double, double> get_wvel(const size_t ii) const {
    return std::make_pair(wvel.data()[2 * ii], wvel.data()[2 * ii + 1]);
  }

  std::pair<double, double> get_uvel(const size_t ii) const {
    return std::make_pair(uvel.data()[2 * ii], uvel.data()[2 * ii + 1]);
  }

  std::pair<double, double> get_vvel(const size_t ii) const {
    return std::make_pair(vvel.data()[2 * ii], vvel.data()[2 * ii + 1]);
  }

  void set_press(const size_t ii, const double p) {
    auto r = press.mutable_unchecked<1>();
    r(ii) = p;
  }

  void set_temp(const size_t ii, const double t) {
    auto r = temp.mutable_unchecked<1>();
    r(ii) = t;
  }

  void set_qvap(const size_t ii, const double qv) {
    auto r = qvap.mutable_unchecked<1>();
    r(ii) = qv;
  }

  void set_qcond(const size_t ii, const double qc) {
    auto r = qcond.mutable_unchecked<1>();
    r(ii) = qc;
  }
};

#endif  // LIBS_CLEO_PYTHON_BINDINGS_COUPLDYN_NUMPY_NUMPY_DYNAMICS_HPP_
