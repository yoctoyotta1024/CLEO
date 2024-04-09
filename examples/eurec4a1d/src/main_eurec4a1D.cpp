/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_eurec4a1D.cpp
 * Project: src
 * Created Date: Tuesday 9th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 9th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * runs the CLEO super-droplet model (SDM) for eurec4a 1-D rainshaft example.
 * After make/compiling, execute for example via:
 * ./src/eurec4a1D ../src/config/config.txt
 */

#include <Kokkos_Core.hpp>
#include <array>
#include <cmath>
#include <concepts>
#include <iostream>
#include <stdexcept>
#include <string_view>
