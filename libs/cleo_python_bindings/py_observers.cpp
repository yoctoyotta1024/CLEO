/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_observers.cpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
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

pyobserver::gridboxes create_gridboxes_observer(const unsigned int interval,
                                                SimpleDataset<FSStore> &dataset,
                                                const size_t maxchunk, const size_t ngbxs);

pyobserver::superdrops create_superdrops_observer(const unsigned int interval,
                                                  SimpleDataset<FSStore> &dataset, FSStore &store,
                                                  const size_t maxchunk);

void pyNullObserver(py::module &m) {
  py::class_<pyca::obs_null>(m, "NullObserver")
      .def(py::init())
      .def("next_obs", &pyca::obs_null::next_obs, py::arg("t_mdl"));
}

void pyObserver(py::module &m) {
  py::class_<pyobserver::obs>(m, "Observer")
      .def(py::init<pyobserver::obs0123456, pyobserver::obs7, pyobserver::mo01234567>())
      .def("next_obs", &pyobserver::obs::next_obs, py::arg("t_mdl"));
}

void pycreate_observer(py::module &m) {
  m.def("pycreate_observer", &create_observer,
        "returns type of Observer suitable for KiD test case", py::arg("config"), py::arg("tsteps"),
        py::arg("dataset"), py::arg("store"));
}

pyobserver::obs create_observer(const Config &config, const Timesteps &tsteps,
                                SimpleDataset<FSStore> &dataset, FSStore &store) {
  const auto enable_observers = config.get_python_bindings().enable_observers;
  const auto obsstep = tsteps.get_obsstep();
  const auto maxchunk = config.get_maxchunk();
  const auto ngbxs = config.get_ngbxs();

  if (!enable_observers.gbxindex) {
    throw std::invalid_argument("gbxindex observer cannot be turned off");
  }
  const Observer auto obs0 = GbxindexObserver(dataset, store, maxchunk, ngbxs);

  auto time_interval = LIMITVALUES::uintmax;
  if (enable_observers.time) {
    time_interval = obsstep;
  }
  const Observer auto obs1 =
      TimeObserver(time_interval, dataset, store, maxchunk, &step2dimlesstime);

  auto totnsupers_interval = LIMITVALUES::uintmax;
  if (enable_observers.totnsupers) {
    totnsupers_interval = obsstep;
  }
  const Observer auto obs2 = TotNsupersObserver(totnsupers_interval, dataset, store, maxchunk);

  auto massmoms_interval = LIMITVALUES::uintmax;
  if (enable_observers.massmoms) {
    massmoms_interval = obsstep;
  }
  const Observer auto obs3 =
      MassMomentsObserver(massmoms_interval, dataset, store, maxchunk, ngbxs);

  auto rainmassmoms_interval = LIMITVALUES::uintmax;
  if (enable_observers.rainmassmoms) {
    rainmassmoms_interval = obsstep;
  }
  const Observer auto obs4 =
      MassMomentsRaindropsObserver(rainmassmoms_interval, dataset, store, maxchunk, ngbxs);

  auto gridboxes_interval = LIMITVALUES::uintmax;
  if (enable_observers.gridboxes) {
    gridboxes_interval = obsstep;
  }
  const Observer auto obs5 =
      create_gridboxes_observer(gridboxes_interval, dataset, maxchunk, ngbxs);

  auto superdrops_interval = LIMITVALUES::uintmax;
  if (enable_observers.superdrops) {
    superdrops_interval = obsstep;
  }
  const Observer auto obs6 =
      create_superdrops_observer(superdrops_interval, dataset, store, maxchunk);

  auto precip_interval = LIMITVALUES::uintmax;
  if (enable_observers.precip) {
    precip_interval = obsstep;
  }
  const Observer auto obs7 =
      MonitorPrecipitationObserver(precip_interval, dataset, store, maxchunk, ngbxs);

  return obs0 >> obs1 >> obs2 >> obs3 >> obs4 >> obs5 >> obs6 >> obs7;
}

pyobserver::gridboxes create_gridboxes_observer(const unsigned int interval,
                                                SimpleDataset<FSStore> &dataset,
                                                const size_t maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<SimpleDataset<FSStore>> auto thermo =
      CollectThermo(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<SimpleDataset<FSStore>> auto windvel =
      CollectWindVel(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<SimpleDataset<FSStore>> auto nsupers =
      CollectNsupers(dataset, maxchunk, ngbxs);

  const CollectDataForDataset<SimpleDataset<FSStore>> auto collect_gbxdata =
      nsupers >> windvel >> thermo;
  return WriteToDatasetObserver(interval, dataset, collect_gbxdata);
}

pyobserver::superdrops create_superdrops_observer(const unsigned int interval,
                                                  SimpleDataset<FSStore> &dataset, FSStore &store,
                                                  const size_t maxchunk) {
  CollectDataForDataset<SimpleDataset<FSStore>> auto sdid = CollectSdId(dataset, maxchunk);
  CollectDataForDataset<SimpleDataset<FSStore>> auto sdgbxindex =
      CollectSdgbxindex(dataset, maxchunk);
  CollectDataForDataset<SimpleDataset<FSStore>> auto xi = CollectXi(dataset, maxchunk);
  CollectDataForDataset<SimpleDataset<FSStore>> auto radius = CollectRadius(dataset, maxchunk);
  CollectDataForDataset<SimpleDataset<FSStore>> auto msol = CollectMsol(dataset, maxchunk);
  CollectDataForDataset<SimpleDataset<FSStore>> auto coord3 = CollectCoord3(dataset, maxchunk);
  CollectDataForDataset<SimpleDataset<FSStore>> auto coord1 = CollectCoord1(dataset, maxchunk);
  CollectDataForDataset<SimpleDataset<FSStore>> auto coord2 = CollectCoord2(dataset, maxchunk);

  const auto collect_sddata =
      coord1 >> coord2 >> coord3 >> msol >> radius >> xi >> sdgbxindex >> sdid;
  return SuperdropsObserver(interval, dataset, store, maxchunk, collect_sddata);
}
