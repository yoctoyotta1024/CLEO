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

// void CartesianMaps::set_nghbrsmaps(const unsigned int coord,
//                                    kokkos_uintmap h_back,
//                                    kokkos_pairmap h_forward)
// {
//   if (coord == 3)
//   {
//     to_back_coord3nghbr = h_back;
//     to_forward_coord3nghbr = h_forward;
//   }
//   else if (coord == 1)
//   {
//     to_back_coord1nghbr = h_back;
//     to_forward_coord1nghbr = h_forward;
//   }
//   else if (coord == 2)
//   {
//     to_back_coord2nghbr = h_back;
//     to_forward_coord2nghbr = h_forward;
//   }
// }