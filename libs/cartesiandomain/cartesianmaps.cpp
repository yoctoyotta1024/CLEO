/*
 * ----- CLEO -----
 * File: cartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 1st November 2023
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

void set_ndims(const size_t dim3,
               const size_t dim1,
               const size_t dim2)
/* sets dimensions (ie. number of gridboxes)
in [coord3, coord1, coord2] directions */
{
  auto h_supers = Kokkos::create_mirror_view(ndims); // mirror ndims in case view is on device memory
  Kokkos::deep_copy(h_supers, supers);

  h_ndims.at(0) = dim3;
  h_ndims.at(1) = dim1;
  h_ndims.at(2) = dim2;
}