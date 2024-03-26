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
 * Last Modified: Tuesday 26th March 2024
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
#include <unordered_map>
#include <vector>

#include "./zarr/dataset.hpp"
#include "./zarr/fsstore.hpp"
// #include "./zarr/zarr_array.hpp"

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

void test_zarr_array(FSStore &store);

template <typename Store>
void test_dataset(Dataset<Store> &dataset);

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    const std::filesystem::path basedir("/home/m/m300950/CLEO/roughpaper/build/bin/dataset.zarr");

    // auto store = FSStore(basedir);
    // test_zarr_array(store);

    auto store = FSStore(basedir);
    auto dataset = Dataset(store);
    test_dataset(dataset);
  }
  Kokkos::finalize();
}

void test_1dzarr(FSStore &store, const viewh_type data, const std::string_view name,
                 const std::vector<size_t> &chunkshape) {
  // create array
  const auto dtype = std::string_view("<f8");
  const auto is_backend = false;
  auto zarr = ZarrArray<FSStore, double>(store, name, dtype, chunkshape, is_backend);

  // output data to array
  zarr.write_to_zarr_array(data);
}

void test_multidzarr(FSStore &store, const viewh_type data, const std::string_view name,
                     const std::vector<size_t> &chunkshape,
                     const std::vector<size_t> &reduced_arrayshape) {
  // create array
  const auto dtype = std::string_view("<f8");
  const auto is_backend = false;
  auto zarr =
      ZarrArray<FSStore, double>(store, name, dtype, chunkshape, is_backend, reduced_arrayshape);

  // output data to array
  zarr.write_to_zarr_array(data);
}

void test_zarr_array(FSStore &store) {
  /* arrays of data returned by observer (maybe on device) */
  auto h_data = observer();

  /* 1 - D test */
  test_1dzarr(store, h_data, "r1d_8", std::vector<size_t>({8}));    // shape = [8], chunks = 0
  test_1dzarr(store, h_data, "r1d_2", std::vector<size_t>({2}));    // shape = [8], chunks = 0,1,2,3
  test_1dzarr(store, h_data, "r1d_6", std::vector<size_t>({6}));    // shape = [8], chunks = 0,1
  test_1dzarr(store, h_data, "r1d_11", std::vector<size_t>({11}));  // shape = [8], chunks = 0

  /* 2-D test */
  test_multidzarr(store, h_data, "m2d_4p2", std::vector<size_t>({4, 2}),
                  std::vector<size_t>({2}));  // shape = [4,2], chunks = 0.0
  test_multidzarr(store, h_data, "m2d_2p2", std::vector<size_t>({2, 2}),
                  std::vector<size_t>({2}));  // shape = [4,2], chunks = 0.0, 1.0
  test_multidzarr(store, h_data, "m2d_4p1", std::vector<size_t>({4, 1}),
                  std::vector<size_t>({2}));  // shape = [4,2], chunks = 0.0, 0.1
  test_multidzarr(store, h_data, "m2d_8p1", std::vector<size_t>({8, 1}),
                  std::vector<size_t>({2}));  // shape = [8,1], chunks = 0.0 WARNING
  test_multidzarr(store, h_data, "m2d_5p2", std::vector<size_t>({5, 2}),
                  std::vector<size_t>({2}));  // shape = [5,2], chunks = 0.0 WARNING

  test_multidzarr(store, h_data, "m2d_3p1", std::vector<size_t>({3, 1}),
                  std::vector<size_t>({2}));  // shape = [5,2], chunks = 0.0, 0.1, 1.0 WARNING
  test_multidzarr(store, h_data, "m2d_11p1", std::vector<size_t>({11, 1}),
                  std::vector<size_t>({2}));  // shape = [8,1], chunks = 0.0 WARNING
  test_multidzarr(store, h_data, "m2d_3p2", std::vector<size_t>({3, 2}),
                  std::vector<size_t>({2}));  // shape = [5,1], chunks = 0.0, 1.0 WARNING

  test_multidzarr(store, h_data, "n2d_5p1", std::vector<size_t>({5, 1}),
                  std::vector<size_t>({1}));  // shape = [8,1], chunks = 0.0, 1.0
  test_multidzarr(store, h_data, "n2d_8p1", std::vector<size_t>({8, 1}),
                  std::vector<size_t>({1}));  // shape = [8,1], chunks = 0.0
  test_multidzarr(store, h_data, "n2d_11p1", std::vector<size_t>({11, 1}),
                  std::vector<size_t>({1}));  // shape = [8,1], chunks = 0.0

  test_multidzarr(store, h_data, "n3d_4p2p1", std::vector<size_t>({4, 2, 1}),
                  std::vector<size_t>({2, 1}));  // shape = [4,2,1], chunks = 0.0.0
  test_multidzarr(store, h_data, "n3d_2p2p1", std::vector<size_t>({2, 2, 1}),
                  std::vector<size_t>({2, 1}));  // shape = [4,2,1], chunks = 0.0.0, 1.0.0
  test_multidzarr(store, h_data, "n3d_1p1p1", std::vector<size_t>({1, 1, 1}),
                  std::vector<size_t>({2, 2}));  // shape = [2,2,2], chunks = (many)

  /* 3-D Test */
  test_multidzarr(store, h_data, "n3d_1p3p2", std::vector<size_t>({1, 3, 2}),
                  std::vector<size_t>({3, 2}));  // shape = [2,3,2], chunks = 1.0.0, 1.0.0 WARNING
  test_multidzarr(store, h_data, "n3d_1p3p1", std::vector<size_t>({1, 3, 1}),
                  std::vector<size_t>({3, 1}));  // shape = [3,3,1] chunks = 0.0.0 -> 2.0.0 WARNING
  test_multidzarr(store, h_data, "n3d_2p3p1", std::vector<size_t>({2, 3, 1}),
                  std::vector<size_t>({3, 1}));  // shape = [4,3,1] chunks = 0.0.0, 1.0.0 WARNING
}

XarrayZarrArray<FSStore, double> make_xarrayzarr(
    FSStore &store, const std::unordered_map<std::string, size_t> &datasetdims,
    const std::string_view name, const std::string_view units, const double scale_factor,
    const std::vector<size_t> &chunkshape, const std::vector<std::string> &dimnames) {
  const auto dtype = std::string_view("<f8");
  return XarrayZarrArray<FSStore, double>(store, datasetdims, name, units, dtype, scale_factor,
                                          chunkshape, dimnames);
}

template <typename Store>
void test_dataset(Dataset<Store> &dataset) {
  // arrays of h_data returned by observer (maybe on device)
  auto h_data = observer();

  dataset.add_dimension({"SdId", 0});
  auto xzarr = dataset.template create_array<double>("xarray_r1d", "m", "<f8", 1e-6, {6}, {"SdId"});

  dataset.set_dimension({"SdId", 10});
  dataset.write_to_array(xzarr, h_data);

  dataset.set_dimension({"SdId", 8});
  dataset.write_arrayshape(xzarr);
}
