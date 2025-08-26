/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: py_observers.hpp
 * Project: cleo_python_bindings
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Python bindings to various different CLEO's Observers instantiations
 */

#ifndef LIBS_CLEO_PYTHON_BINDINGS_PY_OBSERVERS_HPP_
#define LIBS_CLEO_PYTHON_BINDINGS_PY_OBSERVERS_HPP_

#include <pybind11/pybind11.h>

#include <stdexcept>

#include "../cleoconstants.hpp"
#include "./cleo_python_bindings_aliases.hpp"
#include "configuration/config.hpp"
#include "initialise/timesteps.hpp"
#include "observers/collect_data_for_simple_dataset.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/gbxindex_observer.hpp"
#include "observers/massmoments_observer.hpp"
#include "observers/nsupers_observer.hpp"
#include "observers/observers.hpp"
#include "observers/sdmmonitor/do_sdmmonitor_obs.hpp"
#include "observers/sdmmonitor/monitor_precipitation_observer.hpp"
#include "observers/state_observer.hpp"
#include "observers/superdrops_observer.hpp"
#include "observers/time_observer.hpp"
#include "observers/totnsupers_observer.hpp"
#include "zarr/fsstore.hpp"
#include "zarr/simple_dataset.hpp"

namespace py = pybind11;
namespace pyca = cleo_python_bindings_aliases;

void pyNullObserver(py::module &m);
void pyObserver(py::module &m);
void pycreate_observer(py::module &m);

#endif  // LIBS_CLEO_PYTHON_BINDINGS_PY_OBSERVERS_HPP_
