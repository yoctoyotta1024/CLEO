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
 * Last Modified: Thursday 26th June 2025
 * Modified By: CB
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

/* assumes MPI was initialised aleady, e.g. in python with ```from mpi4py import MPI``` */
void pycleo_initialize(const Config &config) {
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  if (comm_size > 1) {
    std::cout << "ERROR: The current example is not prepared"
              << " to be run with more than one MPI process" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (!Kokkos::is_initialized()) {
    Kokkos::initialize(config.get_kokkos_initialization_settings());
    Kokkos::print_configuration(std::cout);
    const auto err = std::atexit(pycleo_finalize);
    if (err) {
      std::cerr << "Failed atexit(pycloe_finalize) in pycleo_initialize()\n";
      return;
    }
  }
}
