/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: wip_fsstore_array.hpp
 * Project: zarr
 * Created Date: Tuesday 12th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 20th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 */

#ifndef ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_
#define ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Pair.hpp>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "./fsstore.hpp"

/**
 * @brief Write attributes string to a store under a .zattrs key.
 *
 * Write some data under .zattrs key in store for an array called 'name'. The key and attrs data
 * could be anything, but for example .zattrs could be a json file in a file system store
 * (see FSStore) for the extra metadata which must exist in order to make xarray and netCDF
 * happy when opening a Zarr dataset, e.g. by naming the dimensions of the
 * "{\"_ARRAY_DIMENSIONS\": [\"dimension_name\"]}";.
 *
 * @tparam Store The type of the store object where the metadata will be written.
 * @param store The store object where the metadata will be written.
 * @param name The name under which the .zarray key will be stored in the store.
 * @param metadata The metadata to write for the .zarray key.
 */
template <typename Store>
inline void write_zattrs_json(Store& store, std::string_view name, std::string_view attrs) {
  store[std::string(name) + "/.zattrs"] = attrs;
}

/* converts vector of strings, e.g. for names of dimensions, into a single list
written as a string */
inline std::string vecstr_to_string(const std::vector<std::string>& dims) {
  auto dims_str = std::string{"["};
  for (const auto& d : dims) {
    dims_str += "\"" + d + "\",";
  }
  dims_str.pop_back();  // delete last ","
  dims_str += "]";
  return dims_str;
}

// TODO(CB) move xarray metadata to dataset
// /* make string of zattrs attribute information for array in zarr store */
// const auto arrayattrs = std::string(
//   "{\n"
//   "  \"_ARRAY_DIMENSIONS\": " +
//   vecstr_to_string(dims) +                // names of each dimension of array
//   ",\n"
//   "  \"units\": " +
//   "\"" + std::string(units) + "\"" +    // units of coordinate being stored
//   ",\n"
//   "  \"scale_factor\": " +
//   std::to_string(scale_factor) +        // scale_factor of data
//   "\n}");

// write_zattrs_json(store, name, arrayattrs);

#endif  // ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_
