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

class ZarrArrayViaBuffer {
 public:
  std::array<double, 1000> buffer;

  ZarrArrayViaBuffer() : buffer{} {};

  ~ZarrArrayViaBuffer() {
    std::cout << "flushing buffer to output\n";
  };

  void write_array(const std::vector<double> &data) {
    std::cout << "writing data to buffer / output\n";
    std::cout << "size: " << data.size() << "\n";
  };
};

#endif    // ROUGHPAPER_ZARR_OUTPUT_HPP_
