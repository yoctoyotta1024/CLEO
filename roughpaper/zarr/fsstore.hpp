/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: fsstore.hpp
 * Project: zarr
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 18th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Class for writing memory in a a file system store under a given key.
 */


#ifndef ROUGHPAPER_ZARR_FSSTORE_HPP_
#define ROUGHPAPER_ZARR_FSSTORE_HPP_

#include <filesystem>
#include <string_view>
#include <span>
#include <fstream>
#include <iostream>

#include "./store_accessor.hpp"

/* A file system (with root in 'basedir' directory) obeying Zarr
version 2 requirements for a Store. Store contins a series
of key, values where values may be data arrays or groups in the store.
data for a given key is written to the store via the functions
in StoreAccessor */
class FSStore {
 private:
  const std::filesystem::path basedir;

 public:
  explicit FSStore(const std::filesystem::path basedir) : basedir(basedir) {}

  StoreAccessor<FSStore> operator[](const std::string_view key) { return { *this, key }; }

  /* write function called by StoreAccessor once data has been
  converted into a vector of unsigned integer types */
  inline bool write(const std::string_view key, const std::span<const uint8_t> buffer);
};

/* write function called by StoreAccessor
once data has been converted into a vector
of unsigned integer types TODO(CB): move to .cpp file*/
inline bool FSStore::write(std::string_view key, std::span<const uint8_t> buffer) {
  auto path = basedir / key;
  auto mode = std::ios::out | std::ios::binary;
  std::ofstream out(path, mode);

  if (!out.good()) {
    std::cout << "couldn't open " << path << ",\n "
      << "making directory " << path.parent_path() << "\n";
    std::filesystem::create_directories(path.parent_path());
    out.open(path, mode);
  }

  if (!out.good()) {
    std::cout << "can't write to " << path << "\n";
    return false;
  }

  out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
  return true;
}

#endif   // ROUGHPAPER_ZARR_FSSTORE_HPP_
