/*
 * ----- CLEO -----
 * File: predcorrmotion.hpp
 * Project: gridboxes
 * Created Date: Tuesday 19th December 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 19th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * -----
 * File Description:
 * Generic struct satisfying Motion concept for 
 * a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * updating gbx according to templated functions
 */


#ifndef PREDCORRMOTION_HPP 
#define PREDCORRMOTION_HPP 

#include <functional>
#include <cassert>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "./predcorr.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"

#endif // PREDCORRMOTION_HPP  