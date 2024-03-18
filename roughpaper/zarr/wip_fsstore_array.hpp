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
 * Last Modified: Monday 18th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 */

#ifndef ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_
#define ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_DualView.hpp>

#include "./fsstore.hpp"

using HostSpace = Kokkos::DefaultHostExecutionSpace;     // TODO(CB) (re-)move definitions
using kkpair_size_t = Kokkos::pair<size_t, size_t>;      // TODO(CB) (re-)move definitions

/* converts vector of strings, e.g. for names of dimensions, into a single list
written as a string */
inline std::string vecstr_to_string(const std::vector<std::string> &dims) {
  auto dims_str = std::string{ "[" };
  for (const auto& d : dims) { dims_str += "\"" + d + "\","; }
  dims_str.pop_back();    // delete last ","
  dims_str += "]";
  return dims_str;
}

#endif    // ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_
