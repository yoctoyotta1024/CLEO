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

#include "zarr2/dataset.hpp"
#include "zarr2/fsstore.hpp"

using viewh_type = Kokkos::View<double *, Kokkos::HostSpace>;  // view of doubles data

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    const std::filesystem::path basedir("/home/m/m300950/CLEO/roughpaper/build/bin/dataset.zarr");
    auto store = FSStore(basedir);
    auto dataset = Dataset(store);

    auto state_observer = StateObserver(obsstep, store, maxchunk, ngbxs);
  }
  Kokkos::finalize();
}

template <typename Store>
void test_dataset(Dataset<Store> &dataset) {
  // arrays of h_data returned by observer (maybe on device)
  auto h_data = observer();

  dataset.add_dimension({"SdId", 0});
  auto xzarr = dataset.template create_array<double>("radius", "m", "<f8", 1e-6, {6},
                                                     {"SdId"});  // shape = [0], chunks = 0,1

  dataset.set_dimension({"SdId", 8});
  dataset.write_to_array(xzarr, h_data);  // shape = [8], chunks = 0,1

  dataset.set_dimension({"SdId", 10});
  dataset.write_arrayshape(xzarr);  // shape = [10], chunks = 0,1
}
