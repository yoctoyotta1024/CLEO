/*
 * ----- CLEO -----
 * File: cartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 2nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions related to creating and using maps to convert
 * between a gridbox indexes and domain coordinates for a 
 * cartesian C grid
 */


#include "./cartesianmaps.hpp"

// TODO 

void CartesianMaps::set_ndims(const size_t ndim3,
                              const size_t ndim1,
                              const size_t ndim2)
{
  // auto h_ndims = Kokkos::create_mirror_view(ndims); // mirror ndims in case view is on device memory

  // h_ndims(0) = ndim3;
  // h_ndims(1) = ndim1;
  // h_ndims(2) = ndim2;

  // Kokkos::deep_copy(ndims, h_ndims);
}

void CartesianMaps::
    set_boundsmaps_via_copy(const kokkos_pairmap::HostMirror h_3,
                            const kokkos_pairmap::HostMirror h_1,
                            const kokkos_pairmap::HostMirror h_2)
{
  // Kokkos::deep_copy(to_coord3bounds, h_3);
  // Kokkos::deep_copy(to_coord1bounds, h_1);
  // Kokkos::deep_copy(to_coord2bounds, h_2);
}

void CartesianMaps::
    set_nghbrsmaps_via_copy(const unsigned int coord,
                            const kokkos_uintmap::HostMirror h_back,
                            const kokkos_pairmap::HostMirror h_forward)
{
  // if (coord == 3)
  // {
  //   Kokkos::deep_copy(to_back_coord3nghbr, h_back);
  //   Kokkos::deep_copy(to_forward_coord3nghbr, h_forward);
  // }
  // else if (coord == 1)
  // {
  //   Kokkos::deep_copy(to_back_coord1nghbr, h_back);
  //   Kokkos::deep_copy(to_forward_coord1nghbr, h_forward);
  // }
  // else if (coord == 2)
  // {
  //   Kokkos::deep_copy(to_back_coord2nghbr, h_back);
  //   Kokkos::deep_copy(to_forward_coord2nghbr, h_forward);
  // }
}