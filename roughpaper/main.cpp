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
 * Last Modified: Wednesday 20th March 2024
 * Last Modified: Wednesday 20th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * rough paper for checking small things
 */

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <iostream>
#include <vector>

#include "./zarr/fsstore.hpp"
#include "./zarr/zarr_array.hpp"

using viewh_type = Kokkos::View<double *, Kokkos::HostSpace>;  // view of doubles data

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

    // auto zarr = ZarrArray<FSStore, double>(store, "radius", "<f8", std::vector<size_t>({9}));
    // auto zarr = ZarrArray<FSStore, double>(store, "massmom", "<f8", std::vector<size_t>({3, 1}),
    //  std::vector<size_t>({2}));
    auto zarr = ZarrArray<FSStore, double>(store, "test3d", "<f8", std::vector<size_t>({1, 4, 1}),
                                           std::vector<size_t>({8, 1}));

    // arrays of data returned by observer (maybe on device)
    auto data = observer();

    // output data to zarr arrays via buffer
    zarr.write_to_array(data);
  }
  Kokkos::finalize();
}
