/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: index_test.cpp
 * Project: roughpaper
 * Created Date: Friday 20th September 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 */

#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  int grid_points = std::stoi(argv[1]);
  std::cout << "grid_points: " << grid_points << "\n";

  size_t vertical_levels = 1;
  size_t ndims_north = 3;
  size_t ndims_east = 3;

  switch (grid_points) {
    case 0:
      /* handles fields defined on centres of cells of horizontal
      (2-D) grid (grid_points = CENTRES), e.g. temp or wvel
       */
      for (size_t j = 0; j < ndims_north; j++) {
        for (size_t i = 0; i < ndims_east; i++) {
          for (size_t k = 0; k < vertical_levels; k++) {
            auto vertical_idx = k;
            auto source_idx = ndims_east * j + i;
            auto ii = (ndims_east * j + i) * vertical_levels + k;
            std::cout << i << ", " << j << " -> " << source_idx << "\n";
          }
        }
      }
      return 0;
    case 1:
      /* handles fields defined on longitude edges of cells of
      horizontal grid (grid_points = LONGITUDE_EDGES), e.g uvel
      */
      ++ndims_east;

      for (size_t j = 0; j < ndims_north; j++) {
        for (size_t i = 0; i < ndims_east; i++) {
          for (size_t k = 0; k < vertical_levels; k++) {
            auto vertical_idx = k;
            auto source_idx = (ndims_east - 1) * j + ndims_east * j;
            source_idx += std::min(2 * i + 1, 2 * ndims_east - 2);
            auto ii = (ndims_east * j + i) * vertical_levels + k;
            std::cout << i << ", " << j << " -> " << source_idx << "\n";
          }
        }
      }
      return 0;

    case 2:
      /* handles fields defined on latitude edges of cells of
      horizontal grid (grid_points = LATITUDE_EDGES), e.g. vvel
      */
      ++ndims_north;

      for (size_t j = 0; j < ndims_north; j++) {
        for (size_t i = 0; i < ndims_east; i++) {
          for (size_t k = 0; k < vertical_levels; k++) {
            auto vertical_idx = k;
            auto source_idx = ndims_east * j + (ndims_east + 1) * j;
            if (j < ndims_north - 1) {
              source_idx += 2 * i;
            } else {
              source_idx += i;
            }
            auto ii = (ndims_east * j + i) * vertical_levels + k;
            std::cout << i << ", " << j << " -> " << source_idx << "\n";
          }
        }
      }
      return 0;
  }
}
