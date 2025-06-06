/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: pycleo.hpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 10th June 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Entry point for creating CLEO's python bindngs module
 */

#ifndef LIBS_PYCLEO_PYCLEO_HPP_
#define LIBS_PYCLEO_PYCLEO_HPP_

#include <mpi.h>
#include <pybind11/pybind11.h>

#include <iostream>

#include "./py_cartesiandomain.hpp"
#include "./py_gridboxes.hpp"
#include "./py_initialise.hpp"
#include "./py_observers.hpp"
#include "./py_runcleo.hpp"
#include "./py_superdrops.hpp"

namespace py = pybind11;

int test_python_bindings(const int i, const int j);
void pycleo_initialize();

void inline pycleo_finalize() {
  Kokkos::finalize();
  // MPI_Finalize();
}

PYBIND11_MODULE(pycleo, m) {
  m.doc() = "Python bindings for selected parts of CLEO's libraries";

  m.def("test_python_bindings", &test_python_bindings, "test function for CLEO example",
        py::arg("i"), py::arg("j"));

  m.def("pycleo_initialize", &pycleo_initialize,
        "necessary initialisation before running CLEO via python ");

  m.def("pycleo_finalize", &pycleo_finalize,
        "necessary finalisation after running CLEO via python ");

  /* initialisation/configuration*/
  pyConfig(m);
  pyTimesteps(m);
  pycreate_timesteps(m);
  pyInitSupersFromBinary(m);
  pyInitSupersFromBinaryParams(m);

  /* superdroplets */
  pySupersInDomain(m);
  pycreate_supers(m);

  /* maps */
  pyCartesianMaps(m);
  pycreate_cartesian_maps(m);

  /* microphyiscs */
  pyNullMicrophysicalProcess(m);

  /* motion */
  pyNullMotion(m);

  /* boundary conditions */
  pyNullBoundaryConditions(m);

  /* transport */
  pyCartesianTransportAcrossDomain(m);

  /* movement */
  pyCartesianNullMoveSupersInDomain(m);

  /* observers */
  pyNullObserver(m);

  /* sdmmethods */
  pyCartesianNullSDMMethods(m);
}

#endif  // LIBS_PYCLEO_PYCLEO_HPP_
