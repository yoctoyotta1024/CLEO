/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: readbinary.cpp
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

#include "initialise/readbinary.hpp"

GblMetadata::GblMetadata(std::ifstream &file) {
  // read 4 unsigned ints are start of binary file
  file.clear();
  file.seekg(0, std::ios::beg);

  std::vector<unsigned int> uints(4, 0);
  binary_into_buffer<unsigned int>(file, uints);

  d0byte = uints.front();
  charbytes = uints.at(1);
  nvars = uints.at(2);
  mbytes_pervar = uints.at(3);

  const unsigned int offset{
      4 * sizeof(unsigned int)};  // offset from start of file to start of metastring
  metastr = read_global_metastring(file, offset);
}

/* read 'gblmbytes' bytes of file and interpret as string of global metadata
to print to terminal. Return current position in file (after reading) */
std::string GblMetadata::read_global_metastring(std::ifstream &file, const int off) const {
  file.seekg(off, std::ios::beg);

  const size_t nchars = charbytes / sizeof(char);
  std::string metastr(nchars, ' ');
  file.read(&metastr[0], nchars);

  std::cout << "----------------- gridfile global metastring -----------------\n"
            << metastr << "\n--------------------------------------------------------------\n";

  return metastr;
}

VarMetadata::VarMetadata(std::ifstream &file, const int off) {
  file.seekg(off, std::ios::beg);

  std::vector<unsigned int> uints(3, 0);
  binary_into_buffer<unsigned int>(file, uints);

  char chars[2];
  file.read(chars, 2 * sizeof(char));

  double dbl;
  file.read(reinterpret_cast<char *>(&dbl), sizeof(double));

  b0 = uints.front();
  bsize = uints.at(1);
  nvar = uints.back();
  vtype = chars[0];
  units = chars[1];
  scale_factor = dbl;
}

/* open binary file for reading or raise error */
std::ifstream open_binary(const std::filesystem::path filename) {
  std::string filestr = filename.string();
  std::cout << "opening binary file: " << filestr << '\n';
  std::ifstream file(filestr, std::ios::in | std::ios::binary);

  if (!file.is_open()) {
    throw std::invalid_argument("Cannot open " + filestr);
  }

  return file;
}

/* Given a binary file that follows the correct layout,
read and print the global metadata string at the start of the file,
then return a vector containing the metadata that is specific to
each of the variables in the file */
std::vector<VarMetadata> metadata_from_binary(std::ifstream &file) {
  const auto gblmeta = GblMetadata(file);

  unsigned int pos = 4 * sizeof(unsigned int) +
                     gblmeta.charbytes;  // position of 1st byte of variable specific metadata

  std::vector<VarMetadata> mdata(0);
  for (unsigned int i = 0; i < gblmeta.nvars; ++i) {
    mdata.push_back(VarMetadata(file, pos));
    pos += gblmeta.mbytes_pervar;
  }

  return mdata;
}

/* raise error if values in vector 'sizes' are not the same. Useful
for checking if vectors are the same size e.g. for vectors of
SD attributes created from reading initSDsfile and used to
make InitSdsData object */
void check_vectorsizes(const std::vector<size_t> &sizes) {
  const size_t sz0 = sizes.front();
  for (auto sz : sizes) {
    if (sz != sz0) {
      const std::string err("values in 'sizes' vector are not identical");
      throw std::invalid_argument(err);
    }
  }
}
