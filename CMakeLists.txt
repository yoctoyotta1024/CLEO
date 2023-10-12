/*
 * ----- CLEO -----
 * File: CMakeLists.txt
 * Project: CLEOfire
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
 * Copyright (c) 2023 Clara Bayley
 */


cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) # if using Kokkos with NVC++ compiler 

project(project
  LANGUAGES CXX
  DESCRIPTION "CLEO by Clara Bayley and other developers"
  VERSION 4.0
)

add_subdirectory(libs)
add_subdirectory(src)