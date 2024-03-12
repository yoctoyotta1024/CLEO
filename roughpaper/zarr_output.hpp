/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: zarr_output.hpp
 * Project: roughpaper
 * Created Date: Tuesday 12th March 2024
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
 */

#ifndef ROUGHPAPER_ZARR_OUTPUT_HPP_
#define ROUGHPAPER_ZARR_OUTPUT_HPP_

#include <iostream>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

using dualview_type = Kokkos::DualView<double *>;             // dual view of doubles

struct Buffer{
 private:
  size_t fill;
  std::array<double, 10> array;

 public:
  void write_to_buffer(const dualview_type data) {
    std::cout << "buffer fill: " << fill << "\n";
    std::cout << "buffer max: " << array.size() << "\n";
    std::cout << "data to add: " << data.extent(0) << "\n";
  }
};

class ZarrArrayViaBuffer {
 public:
  Buffer buffer;

  ZarrArrayViaBuffer() : buffer{} {};

  ~ZarrArrayViaBuffer() {
    std::cout << "flushing buffer to output\n";
  };

  void write_array(const dualview_type data) {
    std::cout << "writing data to buffer / output\n";
    buffer.write_to_buffer(data);
  };
};

#endif    // ROUGHPAPER_ZARR_OUTPUT_HPP_
