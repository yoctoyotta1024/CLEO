/*
 * ----- CLEO -----
 * File: condensation.cpp
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
 * Functionality related to modelling condensation
 * microphysical process in SDM
 */


#include "./condensation.hpp"

KOKKOS_FUNCTION
void DoCondensation::do_condensation(const unsigned int subt) const
/* enact condensation / evaporation microphysical process */
{
  // std::cout << "cond microphys @@ t = " << subt << "\n";
}