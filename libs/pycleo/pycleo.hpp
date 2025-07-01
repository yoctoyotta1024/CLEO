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
 * Last Modified: Tuesday 1st July 2025
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

#include <cstdlib>
#include <iostream>

#include "./py_cartesiandomain.hpp"
#include "./py_configuration.hpp"
#include "./py_gridboxes.hpp"
#include "./py_initialise.hpp"
#include "./py_observers.hpp"
#include "./py_runcleo.hpp"
#include "./py_superdrops.hpp"
#include "configuration/config.hpp"
#include "coupldyn_numpy/pycoupldyn_numpy.hpp"

namespace py = pybind11;

int test_pycleo(const int i, const int j);

void pycleo_initialize(const Config &config);

void inline pycleo_finalize() {
  Kokkos::finalize();
  // MPI_Finalize();
}

PYBIND11_MODULE(pycleo, m) {
  m.doc() = "Python bindings for selected parts of CLEO's libraries";

  m.def("test_pycleo", &test_pycleo, "test function for CLEO example", py::arg("i"), py::arg("j"));

  m.def("pycleo_initialize", &pycleo_initialize,
        "necessary initialisation before running CLEO via python ");

  m.def("pycleo_finalize", &pycleo_finalize,
        "necessary finalisation after running CLEO via python ");

  /* coupldyn_numpy submodule */
  include_coupldyn_numpy_submodule(m);

  /* initialisation/configuration*/
  pyConfig(m);
  pyTimesteps(m);
  pycreate_timesteps(m);
  pyrealtime2step(m);
  pyInitSupersFromBinary(m);
  pyInitSupersFromBinaryParams(m);
  pyInitGbxsNull(m);

  /* superdroplets */
  pySupersInDomain(m);
  pycreate_supers_from_binary(m);

  /* Gridboxes */
  pycreate_gbxs_cartesian_null(m);
  pyGridboxesDualView(m);

  /* maps */
  pyCartesianMaps(m);
  pycreate_cartesian_maps(m);

  /* microphyiscs */
  pyNullMicrophysicalProcess(m);
  pyAllMicrophysicalProcess(m);
  pycreate_microphysical_process(m);

  /* motion */
  pyNullMotion(m);
  pyCartesianPredCorrMotion(m);
  pycreate_cartesian_predcorr_motion(m);

  /* boundary conditions */
  pyNullBoundaryConditions(m);

  /* transport */
  pyCartesianTransportAcrossDomain(m);

  /* movement */
  pyCartesianNullMoveSupersInDomain(m);
  pyCartesianMoveSupersInDomain(m);

  /* observers */
  pyNullObserver(m);

  /* sdmmethods */
  pyCartesianNullSDMMethods(m);
  pyCartesianSDMMethods(m);
}

#endif  // LIBS_PYCLEO_PYCLEO_HPP_
