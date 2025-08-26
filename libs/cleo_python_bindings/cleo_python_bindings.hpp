/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cleo_python_bindings.hpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Entry point for creating CLEO's python bindngs module
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_CLEO_PYTHON_BINDINGS_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_CLEO_PYTHON_BINDINGS_HPP_

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
#include "./py_zarr.hpp"
#include "configuration/communicator.hpp"
#include "configuration/config.hpp"
#include "coupldyn_numpy/pycoupldyn_numpy.hpp"

namespace py = pybind11;

int test_cleo_python_bindings(const int i, const int j);

void cleo_initialize(const Config &config);

void inline cleo_finalize() { Kokkos::finalize(); }

PYBIND11_MODULE(cleo_python_bindings, m) {
  m.doc() = "Python bindings for selected parts of CLEO's libraries";

  m.def("test_cleo_python_bindings", &test_cleo_python_bindings, "test function for CLEO example",
        py::arg("i"), py::arg("j"));

  m.def("cleo_initialize", &cleo_initialize,
        "necessary initialisation before running CLEO via python", py::arg("config"));

  m.def("cleo_finalize", &cleo_finalize, "necessary finalisation after running CLEO via python ");

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
  pyFSStore(m);
  pySimpleDataset(m);
  pyNullObserver(m);
  pyObserver(m);
  pycreate_observer(m);

  /* sdmmethods */
  pyCartesianNullSDMMethods(m);
  pyCartesianSDMMethods(m);
}

#endif  // LIBS_CLEO_PYTHON_BINDINGS_CLEO_PYTHON_BINDINGS_HPP_
