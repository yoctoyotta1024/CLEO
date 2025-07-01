/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_superdrops.cpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 1st July 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for creating python bindings to various parts of CLEO's superdrops library
 */

#include "./py_superdrops.hpp"

CombinedMicrophysicalProcess<NullMicrophysicalProcess, ConstTstepMicrophysics<DoCondensation>>
create_microphysical_process(const Config &config, const Timesteps &tsteps);

void pyNullMicrophysicalProcess(py::module &m) {
  py::class_<pyca::micro_null>(m, "NullMicrophysicalProcess").def(py::init());
}

void pyAllMicrophysicalProcess(py::module &m) {
  py::class_<pyca::micro_all>(m, "AllMicrophysicalProcess")
      .def(py::init<NullMicrophysicalProcess, ConstTstepMicrophysics<DoCondensation>>());
}

void pycreate_microphysical_process(py::module &m) {
  m.def("pycreate_microphysical_process", &create_microphysical_process,
        "returns type of Microphysical Process", py::arg("config"), py::arg("timesteps"));
}

void pyNullMotion(py::module &m) { py::class_<pyca::mo_null>(m, "NullMotion").def(py::init()); }

CombinedMicrophysicalProcess<NullMicrophysicalProcess, ConstTstepMicrophysics<DoCondensation>>
create_microphysical_process(const Config &config, const Timesteps &tsteps) {
  /* Returns combined microphysical process which behaves like a null process unless
  settings for other processes are defined in config.

  Condensation/evaporation created by default with settings such that it's on_step function never
  returns true. However if the paramaters for the condensation configuration struct are set
  (i.e. maxniters is not a NaNVals::sizet()), the an actual active condensation/evaporation
  process is initialised according to this configuration.
  */
  MicrophysicalProcess auto microphys = NullMicrophysicalProcess{};
  std::cout << "Null microphysical process initialised\n";

  const MicrophysicsFunc auto no_cond = DoCondensation(false, 0.0, 0, 0.0, 0.0, 0.0);
  MicrophysicalProcess auto cond = ConstTstepMicrophysics(LIMITVALUES::uintmax, no_cond);

  const auto c = config.get_condensation();
  if (c.maxniters != NaNVals::sizet()) {
    std::cout << "Adding condensation/evaporation to microphysical process\n";
    cond = Condensation(tsteps.get_condstep(), &step2dimlesstime, c.do_alter_thermo, c.maxniters,
                        c.rtol, c.atol, c.MINSUBTSTEP, &realtime2dimless);
  }

  std::cout << "microphysical processes combined\n";
  return microphys >> cond;
}
