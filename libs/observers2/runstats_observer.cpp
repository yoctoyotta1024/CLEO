/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: runstats_observer.cpp
 * Project: observers2
 * Created Date: Monday 8th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 8th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for making and ouputting
 * statistics related to runtime performance
 * e.g. of timestepping CLEO
 */

#include "observers2/runstats_observer.hpp"

void RunStatsObserver::print_summary() const {
  const auto t_start = stats->t_start;
  const auto t_end = stats->t_end;
  const auto t_stepping = t_end - t_start;
  std::cout << std::fixed << std::setprecision(4) << "\n----- CLEO run complete -----\n"
            << "  Initialisation: " << t_start << "s \n"
            << "  Timestepping: " << t_stepping << "s \n"
            << "  Total run duration: " << t_end << "s \n"
            << "-----------------------------\n";
}

/* open and write some the runtime
statistics to a file called 'statsfile' */
void RunStatsObserver::write_to_file() const {
  std::ofstream file(stats_filename);

  if (file.is_open()) {
    const std::string header(
        "### colums are: name duration/s\n"
        "### ---------------------------\n");

    const auto t_start = stats->t_start;
    const auto t_end = stats->t_end;
    const auto t_stepping = t_end - t_start;
    file << header << "init  " << t_start << "\n"
         << "tstep " << t_stepping << "\n"
         << "total " << t_end << "\n";

    file.close();
  } else {
    throw std::runtime_error("unable to open statsfile");
  }
}
