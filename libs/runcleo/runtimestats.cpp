/*
 * ----- CLEO -----
 * File: runtimestats.cpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for making and ouputting statistics
 * related to runtime whilst timestepping CLEO
 */


#include "./runtimestats.hpp"

void RunStats::summary()
{
  std::cout << "\n ----- CLEO run complete ----- \n"
            << "       Duration: " << time2 << "s \n"
            << "       Initialisation: " << time1 << "s \n"
            << "       Timestepping: " << time2 - time1 << "s \n"
            << "------------------------------- \n";
}