/*
 * ----- CLEO -----
 * File: collisions.cpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for modelling collision
 * microphysical processes in SDM
 * e.g. collision-coalescence
 */

#include "./collisions.hpp"

KOKKOS_FUNCTION
void DoCollisions::do_collisions(const unsigned int subt) const
{
  // std::cout << "coll microphys @@ t = " << subt << "\n";
}