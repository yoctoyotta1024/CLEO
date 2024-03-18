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
 * Last Modified: Monday 18th March 2024
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

#include "./zarr/fsstore_array.hpp"

using viewh_type = Kokkos::View<double *, HostSpace::memory_space>;  // view of doubles data

viewh_type observer() {
  auto h_data = viewh_type("data", 8);

  // initialise data in host view
  h_data(0) = 1.1;
  h_data(1) = 2.2;
  h_data(2) = 3.3;
  h_data(3) = 4.4;
  h_data(4) = 5.5;
  h_data(5) = 6.6;
  h_data(6) = 7.7;
  h_data(7) = 8.8;

  return h_data;
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    const std::filesystem::path basedir("/home/m/m300950/CLEO/roughpaper/build/bin/dataset.zarr");
    auto store = FSStore(basedir);

    // auto zarr = FSStoreArrayViaBuffer<double>(store, std::vector<size_t>({6}), "radius",
    //   "micro-m", 10.0, "<f8", std::vector<std::string>({"sdId"}));

    // auto zarr = FSStoreArrayViaBuffer<double>(store, std::vector<size_t>({5, 1}), "massmom",
    //   "", 1.0, "<f8", std::vector<std::string>({"time", "gbx"}),
    //   std::vector<size_t>({2}));

    auto zarr = FSStoreArrayViaBuffer<double>(store, std::vector<size_t>({12, 2, 1}), "radius3d",
      "", 1.0, "<f8", std::vector<std::string>({"time", "gbx", "sdId"}),
      std::vector<size_t>({2, 1}));

    // arrays of data returned by observer (maybe on device)
    auto data = observer();

    // output data to zarr arrays via buffer
    zarr.write_data_to_zarr_array(data);
  }
  Kokkos::finalize();
}
