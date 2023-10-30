/*
 * ----- CLEO -----
 * File: readbinary.cpp
 * Project: initialise
 * Created Date: Monday 30th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 30th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * tools for reding binary initialisation
 * file e.g. for making gridbox maps or
 * SD initial conditions */


#include "readbinary.hpp"

std::ifstream open_binary(std::string_view filename)
/* open binary file for reading or raise error */
{
  std::string filestr = static_cast<std::string>(filename);
  std::cout << "opening binary file: " << filestr << '\n';
  std::ifstream file(filestr, std::ios::in | std::ios::binary);

  if (!file.is_open())
  {
    throw std::invalid_argument("Cannot open " + filestr);
  }

  return file;
}
