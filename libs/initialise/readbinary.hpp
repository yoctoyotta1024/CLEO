/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: readbinary.hpp
 * Project: initialise
 * Created Date: Monday 30th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * tools for reding binary initialisation
 * file e.g. for making gridbox maps or
 * SD initial conditions */

#ifndef LIBS_INITIALISE_READBINARY_HPP_
#define LIBS_INITIALISE_READBINARY_HPP_

#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

/* Global Metadata is 4 unsigned integers at very
start of binary (used to declare layout of binary file)
as well as string describing data in file */
struct GblMetadata {
  unsigned int d0byte;         // position of first byte of data
  unsigned int charbytes;      // no. bytes of global metadata chars (in string)
  unsigned int nvars;          // no. variables in file
  unsigned int mbytes_pervar;  //  no. bytes of metadata per variable
  std::string metastr;         // description of file contents

  explicit GblMetadata(std::ifstream &file);

  std::string read_global_metastring(std::ifstream &file, const int off) const;
};

/* metadata in file related to a
particular variable (vaR) in that file */
struct VarMetadata {
  unsigned int b0;      // first byte in file containing this var's data
  unsigned int bsize;   // size in bytes of 1 datapoint of this var
  unsigned int nvar;    // no. datapoints of this var
  char vtype;           // char indicating type of this var
  char units;           // char indicating units once data multiplied by scale_factor
  double scale_factor;  // scale factor to re-dimensionalise data

  VarMetadata() {}

  VarMetadata(std::ifstream &file, const int off);
};

/* open binary file for reading or raise error */
std::ifstream open_binary(const std::filesystem::path filename);

/* Given a binary file that follows the correct layout,
read and print the global metadata string at the start of the file,
then return a vector containing the metadata that is specific to
each of the variables in the file */
std::vector<VarMetadata> metadata_from_binary(std::ifstream &file);

/* raise error if values in vector 'sizes' are not the same. Useful
for checking if vectors are the same size e.g. for vectors of
SD attributes created from reading initSDsfile and used to
make InitSdsData object */
void check_vectorsizes(const std::vector<size_t> &sizes);

template <typename T>
void binary_into_buffer(std::ifstream &file, std::vector<T> &buffer) {
  file.read(reinterpret_cast<char *>(buffer.data()), buffer.size() * sizeof(T));
}

/* return vector of data read from ifstream file for
one variable in a binary file given that variable's
metadata is given by the VarMetadata instance, 'varmeta' */
template <typename T>
std::vector<T> vector_from_binary(std::ifstream &file, const VarMetadata &varmeta) {
  file.seekg(varmeta.b0, std::ios::beg);
  std::vector<T> vardata(varmeta.nvar, 0);
  binary_into_buffer<T>(file, vardata);

  return vardata;  // data for variable in binary file given it's metadata
}

#endif  // LIBS_INITIALISE_READBINARY_HPP_
