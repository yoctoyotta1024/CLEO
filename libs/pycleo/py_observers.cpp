/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_observers.cpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 15th July 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for creating python bindings to various different CLEO's Observers instantiations
 */

#include "./py_observers.hpp"

pyobserver::obs create_observer(const Config &config, const Timesteps &tsteps,
                                SimpleDataset<FSStore> &dataset, FSStore &store);

void pyNullObserver(py::module &m) {
  py::class_<pyca::obs_null>(m, "NullObserver")
      .def(py::init())
      .def("next_obs", &pyca::obs_null::next_obs, py::arg("t_mdl"));
}

void pyObserver(py::module &m) {
  py::class_<pyobserver::obs>(m, "Observer")
      .def(py::init<pyobserver::gbx, pyobserver::time, pyobserver::mo>())
      .def("next_obs", &pyobserver::obs::next_obs, py::arg("t_mdl"));
}

void pycreate_observer(py::module &m) {
  m.def("pycreate_observer", &create_observer,
        "returns type of Observer suitable for KiD test case", py::arg("config"), py::arg("tsteps"),
        py::arg("dataset"), py::arg("store"));
}

pyobserver::obs create_observer(const Config &config, const Timesteps &tsteps,
                                SimpleDataset<FSStore> &dataset, FSStore &store) {
  const auto enable_observers = config.get_pycleo().enable_observers;
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();
  const auto ngbxs = config.get_ngbxs();

  auto time_interval = LIMITVALUES::uintmax;
  if (enable_observers.time) {
    time_interval = obsstep;
  }
  if (!enable_observers.gbxindex) {
    throw std::invalid_argument("gbxindex observer cannot be turned off");
  }

  const Observer auto obs1 =
      TimeObserver(time_interval, dataset, store, maxchunk, &step2dimlesstime);
  const Observer auto obs2 = GbxindexObserver(dataset, store, maxchunk, ngbxs);

  return obs2 >> obs1;
}
