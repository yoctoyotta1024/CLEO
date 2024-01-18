/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: longhydroprob.cpp
 * Project: collisionprobs
 * Created Date: Thursday 9th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 22nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for probability of collision-coalescence
 * event between two (real) droplets using the
 * hydrodynamic (i.e. gravitational) kernel
 * according to Simmel et al. 2002's formulation
 * of Long's Hydrodynamic Kernel. Probability
 * calculations are contained in structures
 * that satisfy the requirements of the
 * PairProbability concept (see collisions.hpp)
 */


#include "./longhydroprob.hpp"
