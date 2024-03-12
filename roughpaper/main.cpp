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
#include <vector>

#include <Kokkos_Core.hpp>

#include "./zarr_output.hpp"

int main(int argc, char *argv[]) {
  auto zarr = ZarrArrayViaBuffer();

  Kokkos::initialize(argc, argv);
  {
    // arrays of data returned by observer (maybe on device)
    auto data1 = std::vector<double>(100, 1.1);
    auto data2 = std::vector<double>(1000, 22.22);
    auto data3 = std::vector<double>(5000, 333.333);

    // output data to zarr arrays via buffer
    zarr.write_array(data1);
    zarr.write_array(data2);
    zarr.write_array(data3);
  }
  Kokkos::finalize();
}
