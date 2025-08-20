/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: pycleo.cpp
 * Project: pycleo
 * Created Date: Thursday 5th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for creating CLEO's python bindngs module
 */

#include "./pycleo.hpp"

int test_pycleo(const int i, const int j) {
  std::cout << "Hello World\n";
  return i + j;
}

void pycleo_initialize(const Config &config) {
  /* Initialize Communicator here */
  /* NOTE: call to init_communicator constructor assumes MPI was initialised aleady,
    e.g. in python with ```from mpi4py import MPI``` */
  init_communicator init_comm(0, NULL, config);

  /* Prevent python bindings from running with more than one MPI process */
  const auto comm_size = init_communicator::get_comm_size();
  if (comm_size > 1) {
    throw std::invalid_argument(
        "ERROR: The current example is not prepared to be run with more than one MPI process");
  }

  /* Initialize Kokkos if not already initialised here */
  if (!Kokkos::is_initialized()) {
    Kokkos::initialize(config.get_kokkos_initialization_settings());
    Kokkos::print_configuration(std::cout);
    const auto err = std::atexit(pycleo_finalize);
    if (err) {
      std::cerr << "Failed atexit(pycleo_finalize) in pycleo_initialize()\n";
      return;
    }
  }
}
