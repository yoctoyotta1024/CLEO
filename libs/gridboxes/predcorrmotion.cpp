/*
 * ----- CLEO -----
 * File: predcorrmotion.cpp
 * Project: gridboxes
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 7th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Functionality for motion using predictor-corrector
 * method to update a superdroplet's coordinates given
 * a formula for its terminal velocity and the wind
 * velocity obtained via a simple linear interpolation.
 * Methods follows equations in Grabowski et al. 2018
 */

#include "./predcorrmotion.hpp"

KOKKOS_FUNCTION
void PredCorrMotion::
    update_superdrop_coords(const unsigned int t_sdm) const
{
}