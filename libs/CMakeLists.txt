/*
 * ----- CLEO -----
 * File: CMakeLists.txt
 * Project: libs
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 12th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * CMake for CLEO libraries
 */


cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) # if using Kokkos with NVC++ compiler 

# essential CLEO SDM libaries #
add_subdirectory(initialise)
add_subdirectory(observers)
add_subdirectory(sdmdomain)
add_subdirectory(superdrops)
add_subdirectory(zarrstorage)

# libraries particular to certain CLEO configurations #
add_subdirectory(cvode EXCLUDE_FROM_ALL)
add_subdirectory(thermofromfile EXCLUDE_FROM_ALL)