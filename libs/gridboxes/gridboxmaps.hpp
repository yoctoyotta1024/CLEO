/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: gridboxmaps.hpp
 * Project: gridboxes
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept for maps to convert between a gridbox indexes
 * and domain coordinates for a and type of C grid used
 * by CLEO SDM
 */

#ifndef LIBS_GRIDBOXES_GRIDBOXMAPS_HPP_
#define LIBS_GRIDBOXES_GRIDBOXMAPS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <concepts>

/* concept for GridboxMaps is all types that have
correct signatues for map-like functions */
template <typename GbxMaps>
concept GridboxMaps = requires(GbxMaps gbxmaps, unsigned int idx, size_t s, double& d) {
  { gbxmaps.get_total_global_ngridboxes() } -> std::convertible_to<size_t>;
  { gbxmaps.get_local_ngridboxes() } -> std::convertible_to<size_t>;
  { gbxmaps.get_local_ngridboxes_hostcopy() } -> std::convertible_to<size_t>;

  { gbxmaps.get_gbxarea(idx) } -> std::convertible_to<double>;

  { gbxmaps.get_gbxvolume(idx) } -> std::convertible_to<double>;

  { gbxmaps.coord3bounds(idx) } -> std::convertible_to<Kokkos::pair<double, double>>;
  { gbxmaps.coord1bounds(idx) } -> std::convertible_to<Kokkos::pair<double, double>>;
  { gbxmaps.coord2bounds(idx) } -> std::convertible_to<Kokkos::pair<double, double>>;

  { gbxmaps.coord3backward(idx) } -> std::convertible_to<unsigned int>;
  { gbxmaps.coord3forward(idx) } -> std::convertible_to<unsigned int>;

  { gbxmaps.coord1backward(idx) } -> std::convertible_to<unsigned int>;
  { gbxmaps.coord1forward(idx) } -> std::convertible_to<unsigned int>;

  { gbxmaps.coord2backward(idx) } -> std::convertible_to<unsigned int>;
  { gbxmaps.coord2forward(idx) } -> std::convertible_to<unsigned int>;

  { gbxmaps.global_to_local_gbxindex(s) } -> std::convertible_to<unsigned int>;
  { gbxmaps.local_to_global_gridbox_index(idx) } -> std::convertible_to<size_t>;
  { gbxmaps.get_local_bounding_gridbox_index(idx, d, d, d) } -> std::convertible_to<unsigned int>;
};

#endif  // LIBS_GRIDBOXES_GRIDBOXMAPS_HPP_
