/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_runcleo.cpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 24th June 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for creating python bindings to parts of the runcleo library e.g.
 * to various different CLEO's SDMMethods instantiations
 */

#include "./py_runcleo.hpp"

void pycreate_supers_from_binary(py::module &m) {
  m.def("create_supers_from_binary", &create_supers<InitSupersFromBinary>,
        "returns SupersInDomain instance", py::arg("sdic"), py::arg("gbxindex_max"));
}

void pycreate_gbxs_cartesian_null(py::module &m) {
  m.def("create_gbxs_cartesian_null", &create_gbxs<pyca::map_cart, InitGbxsNull>,
        "returns dualview of Gridboxes instance", py::arg("gbxmaps"), py::arg("gbxic"),
        py::arg("allsupers"));
}

void pyCartesianNullSDMMethods(py::module &m) {
  py::class_<pyca::sdm_cart_null>(m, "CartesianNullSDMMethods")
      .def(py::init<const unsigned int, pyca::map_cart, pyca::micro_null, pyca::move_cart_null,
                    pyca::obs_null>())
      .def_readonly("gbxmaps", &pyca::sdm_cart_null::gbxmaps)
      .def_readonly("obs", &pyca::sdm_cart_null::obs)
      .def("get_couplstep", &pyca::sdm_cart_null::get_couplstep)
      .def("next_couplstep", &pyca::sdm_cart_null::next_couplstep, py::arg("t_mdl"))
      .def(
          "prepare_to_timestep",
          [](const pyca::sdm_cart_null &self, const dualview_gbx gbxs) {
            self.prepare_to_timestep(gbxs.view_device());
          },
          py::arg("gbxs"))
      .def("at_start_step", &pyca::sdm_cart_null::at_start_step, py::arg("t_mdl"), py::arg("gbxs"),
           py::arg("allsupers"))
      .def(
          "run_step",
          [](const pyca::sdm_cart_null &self, const unsigned int t_mdl,
             const unsigned int t_mdl_next, const dualview_gbx gbxs, SupersInDomain &allsupers) {
            self.run_step(t_mdl, t_mdl_next, gbxs.view_device(), allsupers);
          },
          py::arg("t_mdl"), py::arg("t_mdl_next"), py::arg("gbxs"), py::arg("allsupers"));
}

void pyCartesianCondSDMMethods(py::module &m) {
  py::class_<pyca::sdm_cart_cond>(m, "CartesianCondSDMMethods")
      .def(py::init<const unsigned int, pyca::map_cart, pyca::micro_cond, pyca::move_cart_null,
                    pyca::obs_null>())
      .def_readonly("gbxmaps", &pyca::sdm_cart_cond::gbxmaps)
      .def_readonly("obs", &pyca::sdm_cart_cond::obs)
      .def("get_couplstep", &pyca::sdm_cart_cond::get_couplstep)
      .def("next_couplstep", &pyca::sdm_cart_cond::next_couplstep, py::arg("t_mdl"))
      .def(
          "prepare_to_timestep",
          [](const pyca::sdm_cart_cond &self, const dualview_gbx gbxs) {
            self.prepare_to_timestep(gbxs.view_device());
          },
          py::arg("gbxs"))
      .def("at_start_step", &pyca::sdm_cart_cond::at_start_step, py::arg("t_mdl"), py::arg("gbxs"),
           py::arg("allsupers"))
      .def(
          "run_step",
          [](const pyca::sdm_cart_cond &self, const unsigned int t_mdl,
             const unsigned int t_mdl_next, const dualview_gbx gbxs, SupersInDomain &allsupers) {
            self.run_step(t_mdl, t_mdl_next, gbxs.view_device(), allsupers);
          },
          py::arg("t_mdl"), py::arg("t_mdl_next"), py::arg("gbxs"), py::arg("allsupers"));
}
