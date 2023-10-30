/*
 * ----- CLEO -----
 * File: readbinary.hpp
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


#ifndef READBINARY_HPP
#define READBINARY_HPP

// #include <string>
// #include <string_view>
// #include <iostream>
// #include <ios>
// #include <fstream>
// #include <istream>
// #include <vector>
// #include <stdexcept>


std::ifstream open_binary(std::string_view filename);
/* open binary file for reading or raise error */

template <typename T>
std::vector<T> vector_from_binary(std::ifstream &file,
                                  const VarMetadata &varmeta)
/* return vector of data read from ifstream file for
one variable in a binary file given that variable's 
metadata is given by the VarMetadata instance, 'varmeta' */
{
  file.seekg(varmeta.b0, std::ios::beg);
  std::vector<T> vardata(varmeta.nvar, 0);
  binary_into_buffer<T>(file, vardata);

  return vardata; // data for variable in binary file given it's metadata
}


#endif // READBINARY_HPP


