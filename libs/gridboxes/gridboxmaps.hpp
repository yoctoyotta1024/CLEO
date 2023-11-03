/*
 * ----- CLEO -----
 * File: gridboxmaps.hpp
 * Project: gridboxes
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * concept for maps to convert between a gridbox indexes
 * and domain coordinates for a and type of C grid used
 * by CLEO SDM
 */


#ifndef GRIDBOXMAPS_HPP
#define GRIDBOXMAPS_HPP

#include <concepts>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

template <typename GbxMaps>
concept GridboxMaps = requires(GbxMaps gbxmaps,
                               unsigned int idx)
/* concept for GridboxMaps is all types that have
correct signatues for map-like functions */
{ 
  {
    gbxmaps.maps_size()
  } -> std::same_as<size_t>;

  {
    gbxmaps.get_gbxarea(idx)
  } -> std::convertible_to<double>;

  {
    gbxmaps.get_gbxvolume(idx)
  } -> std::convertible_to<double>;

  {
    gbxmaps.coord3bounds(idx)
  } -> std::convertible_to<Kokkos::pair<double, double>>;
  {
    gbxmaps.coord1bounds(idx)
  } -> std::convertible_to<Kokkos::pair<double, double>>;
  {
    gbxmaps.coord2bounds(idx)
  } -> std::convertible_to<Kokkos::pair<double, double>>;

  {
    gbxmaps.coord3backward(idx)
  } -> std::convertible_to<double>;
  {
    gbxmaps.coord3forward(idx)
  } -> std::convertible_to<double>;

  {
    gbxmaps.coord1backward(idx)
  } -> std::convertible_to<double>;
  {
    gbxmaps.coord1forward(idx)
  } -> std::convertible_to<double>;

  {
    gbxmaps.coord2backward(idx)
  } -> std::convertible_to<double>;
  {
    gbxmaps.coord2forward(idx)
  } -> std::convertible_to<double>;
};

#endif // GRIDBOXMAPS_HPP