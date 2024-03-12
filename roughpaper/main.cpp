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
#include <Kokkos_DualView.hpp>

#include "./zarr_output.hpp"

dualview_type observer() {
  auto data = dualview_type("data", 5);

  // initialise data in host view
  auto h_data = data.view_host();
  h_data(0) = 1.1;
  h_data(1) = 2.2;
  h_data(2) = 3.3;
  h_data(3) = 4.4;
  h_data(4) = 5.5;
  data.modify_host();

  // match device data with host
  data.sync_device();

  return data;
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    auto zarr = ZarrArrayViaBuffer(8);

    // arrays of data returned by observer (maybe on device)
    auto data = observer();

    // output data to zarr arrays via buffer
    zarr.write_array(data.view_host());
    std::cout << "--\n";
  }
  Kokkos::finalize();
}
