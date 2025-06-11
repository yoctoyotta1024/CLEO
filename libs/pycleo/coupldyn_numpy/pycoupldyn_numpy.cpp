/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: pycoupldyn_numpy.cpp
 * Project: coupldyn_numpy
 * Created Date: Wednesday 11th June 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 11th June 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for creating CLEO's python bindngs module for the coupldyn_numpy library
 */

#include "./pycoupldyn_numpy.hpp"

int test_python_bindings(const int i, const int j) {
  std::cout << "Hello World\n";
  return i * j;
}
