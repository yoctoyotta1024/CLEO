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
 * Last Modified: Thursday 5th June 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for creating CLEO's python bindngs
 */

#include "./pycleo.hpp"

int test_python_bindings(const int i, const int j) {
  std::cout << "Hello World\n";
  return i + j;
}
