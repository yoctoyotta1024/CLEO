/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: zfsstore.hpp
 * Project: zarr
 * Created Date: Tuesday 12th March 2024
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
 */

#ifndef ROUGHPAPER_ZARR_ZFSSTORE_HPP_
#define ROUGHPAPER_ZARR_ZFSSTORE_HPP_

#include <span>
#include <string>
#include <string_view>
#include <filesystem>
#include <fstream>
#include <iostream>

 /* functions for converting types (e.g. vectors or Kokkos view's of
 unsigned integers or doubles) into vectors of single bytes to
 write to store under a given key. Store can be anything that
 satisfies the zarr storage specifcaiton version 2 */
template <typename Store>
struct StoreAccessor {
  using HostSpace = Kokkos::DefaultHostExecutionSpace;    // TODO(CB) (re-)move definitions

  Store& store;
  std::string_view key;

  /* write range of memory representing uint8_ts to store */
  StoreAccessor& operator=(std::span<const uint8_t> buffer) {
    store.write(key, buffer);
    return *this;
  }

  /* reinterpret range of memory representing string as
  a range of memory representing uint8_ts, then write to store */
  StoreAccessor& operator=(std::string_view buffer) {
    return operator=(
      std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(buffer.data()), buffer.size()));
  }

  /* re-interpret range of memory representing vector of type T as
  a range of memory representing uint8_ts, then write to store */
  template <typename T>
  StoreAccessor& operator=(std::span<const T> buffer) {
    return operator=(std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(buffer.data()),
      buffer.size() * sizeof(T)));
  }

  /* re-interpret range of memory representing vector of type T as
  a range of memory representing uint8_ts, then write to store */
  template <typename T>
  StoreAccessor& operator=(const Kokkos::View<T*, HostSpace::memory_space> buffer) {
    return operator=(std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(buffer.data()),
      buffer.extent(0) * sizeof(T)));
  }
};

/* A file system (with root in 'basedir' directory) obeying Zarr
version 2 requirements for a Store. Store contins a series
of key, values where values may be data arrays or groups in the store.
data for a given key is written to the store via the functions
in StoreAccessor */
class FSStore {
 private:
  const std::filesystem::path basedir;

 public:
  explicit FSStore(const std::filesystem::path basedir) : basedir(basedir) {
    // initialize a zarr group (i.e. dataset)
    const std::string zarr_format("2");   // storage spec. version 2
    const std::string zgroupjson("{\n  \"zarr_format\": " + zarr_format + "\n}");
    (*this)[".zgroup"] = zgroupjson;

    // global metadata (optional)
    (*this)[".zattrs"] =
      "{\n"
      "  \"creator\": \"Clara Bayley\",\n"
      "  \"title\": \"Zarr File System Store for Output Data from CLEO\""
      "\n}";
  }

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

/* write .zarray json file into store for the metadata of an array called 'name'. This metadata MUST
exist in order to decode chunks of an array in a Zarr FSStore. */
inline void write_zarray_json(FSStore& store, std::string_view name,
  std::string_view zarr_metadata) {
  store[std::string(name) + "/.zarray"] = zarr_metadata;
}

/* write .zattr json file into store for the attributes of an array or group called 'name'. attrs
is optional, i.e. not strictly required by Zarr storage specification, but can be useful. For
example attrs can define a group of arrays. Or attrs can define the names of the dimensions of an
array in order to make xarray and netCDF happy e.g. "{\"_ARRAY_DIMENSIONS\": [\"x\"]}"; */
inline void write_zattrs_json(FSStore& store, std::string_view name, std::string_view attrs) {
  store[std::string(name) + "/.zattrs"] = attrs;
}

#endif   // ROUGHPAPER_ZARR_ZFSSTORE_HPP_
