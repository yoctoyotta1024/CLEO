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
 * Last Modified: Thursday 21st March 2024
 * Last Modified: Thursday 21st March 2024
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

void test_1d(FSStore &store, const viewh_type data, const std::string_view name,
             const std::vector<size_t> &chunkshape) {
  // create array
  const auto dtype = std::string_view("<f8");
  auto zarr = ZarrArray<FSStore, double>(store, name, dtype, chunkshape);

  // output data to array
  zarr.write_to_zarr_array(data);
}

void test_multid(FSStore &store, const viewh_type data, const std::string_view name,
                 const std::vector<size_t> &chunkshape,
                 const std::vector<size_t> &reduced_arrayshape) {
  // create array
  const auto dtype = std::string_view("<f8");
  auto zarr = ZarrArray<FSStore, double>(store, name, dtype, chunkshape, reduced_arrayshape);

  // output data to array
  zarr.write_to_zarr_array(data);
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    const std::filesystem::path basedir("/home/m/m300950/CLEO/roughpaper/build/bin/dataset.zarr");
    auto store = FSStore(basedir);

    // arrays of data returned by observer (maybe on device)
    auto data = observer();
    // test_1d(store, data, "r1d_8", std::vector<size_t>({8}));
    // test_1d(store, data, "r1d_2", std::vector<size_t>({2}));
    // test_1d(store, data, "r1d_6", std::vector<size_t>({6}));
    // test_1d(store, data, "r1d_11", std::vector<size_t>({11}));

    // test_multid(store, data, "m2d_4p2", std::vector<size_t>({4, 2}), std::vector<size_t>({2}));
    // test_multid(store, data, "m2d_2p2", std::vector<size_t>({2, 2}), std::vector<size_t>({2}));
    // test_multid(store, data, "m2d_5p2", std::vector<size_t>({5, 2}), std::vector<size_t>({2}));
    // test_multid(store, data, "m2d_8p1", std::vector<size_t>({8, 1}), std::vector<size_t>({2}));
    // test_multid(store, data, "m2d_4p1", std::vector<size_t>({4, 1}), std::vector<size_t>({2}));
    test_multid(store, data, "m2d_3p1", std::vector<size_t>({3, 1}),
                std::vector<size_t>({2}));  // good example of array shape problem
    test_multid(store, data, "m2d_11p1", std::vector<size_t>({11, 1}),
                std::vector<size_t>({2}));  // another example of array shape problem
    test_multid(store, data, "m2d_3p2", std::vector<size_t>({3, 2}),
                std::vector<size_t>({2}));  // another example of array shape problem

    // test_multid(store, data, "n2d_51", std::vector<size_t>({5, 1}), std::vector<size_t>({1}));
    // test_multid(store, data, "n2d_81", std::vector<size_t>({8, 1}), std::vector<size_t>({1}));
  }
  Kokkos::finalize();
}
