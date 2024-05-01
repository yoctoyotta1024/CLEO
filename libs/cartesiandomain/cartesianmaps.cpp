/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 1st May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions related to creating and using maps to convert
 * between a gridbox indexes and domain coordinates for a
 * cartesian C grid
 */

#include "cartesiandomain/cartesianmaps.hpp"

/* on host, throws error if maps are not all
the same size, else returns size of maps */
size_t CartesianMaps::maps_size() const {
  const auto h_ndims(ndims_hostcopy());
  const size_t sz(h_ndims(0) * h_ndims(1) * h_ndims(2) + 1);  // ngbxs + 1 for out of bounds key

  if (to_coord3bounds.size() != sz || to_coord1bounds.size() != sz ||
      to_coord2bounds.size() != sz || to_back_coord3nghbr.size() != sz ||
      to_forward_coord3nghbr.size() != sz || to_back_coord1nghbr.size() != sz ||
      to_forward_coord1nghbr.size() != sz || to_back_coord2nghbr.size() != sz ||
      to_forward_coord2nghbr.size() != sz || to_area.size() != sz || to_volume.size() != sz) {
    throw std::invalid_argument("gridbox maps are not all the same size");
  }

  return sz;
}
