/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Monday 29th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 12th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * rough paper for checking small things
 */

#include <iostream>
#include <array>

#include <Kokkos_Core.hpp>

#include "./zarr_output.hpp"

int main(int argc, char *argv[]) {
  auto zarr = ZarrArrayViaBuffer();

  Kokkos::initialize(argc, argv);
  {
    // arrays of data returned by observer (maybe on device)
    auto data1 = std::array<double, 100>{};
    data1.fill(1.1);

    auto data2 = std::array<double, 1000>{};
    data2.fill(10.0);

    auto data3 = std::array<double, 5000>{};
    data3.fill(222.222);
  }
  Kokkos::finalize();
}
