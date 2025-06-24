/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_runcleo.hpp
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
 * Python bindings to parts of the runcleo library e.g.
 * to various different CLEO's SDMMethods instantiations
 */

#ifndef LIBS_PYCLEO_PY_RUNCLEO_HPP_
#define LIBS_PYCLEO_PY_RUNCLEO_HPP_

#include <pybind11/pybind11.h>

#include "../kokkosaliases.hpp"
#include "./pycleo_aliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_transport_across_domain.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "initialise/init_supers_from_binary.hpp"
#include "initialise/initgbxsnull.hpp"
#include "observers/observers.hpp"
#include "runcleo/creategbxs.hpp"
#include "runcleo/createsupers.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/condensation.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"

namespace py = pybind11;
namespace pyca = pycleo_aliases;

void pycreate_supers_from_binary(py::module &m);

void pycreate_gbxs_cartesian_null(py::module &m);

void pyCartesianNullSDMMethods(py::module &m);
void pyCartesianCondSDMMethods(py::module &m);

#endif  // LIBS_PYCLEO_PY_RUNCLEO_HPP_
