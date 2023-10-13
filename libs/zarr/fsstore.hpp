/*
 * ----- CLEO -----
 * File: fsstore.hpp
 * Project: zarr
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
 */

#ifndef FSSTORE_HPP
#define FSSTORE_HPP

#include <string>
#include <string_view>
#include <span>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>

template <typename Store>
struct StoreAccessor
/* functions for converting types (e.g. vectors of
unsigned integers or doubles) into vectors of single bytes to
write to store under a given key. Store can be anything that
satisfies the zarr storage specifcaiton version 2 */
{
  Store &store;
  std::string_view key;

  StoreAccessor &operator=(std::span<const uint8_t> buffer)
  /* write range of memory representing uint8_ts to store */
  {
    store.write(key, buffer);
    return *this;
  }

  StoreAccessor &operator=(std::string_view buffer)
  /* reinterpret range of memory representing string as
  a range of memory representing uint8_ts, then write to store */
  {
    return operator=(std::span<const uint8_t>(
        reinterpret_cast<const uint8_t *>(buffer.data()),
        buffer.size()));
  }

  template <typename T>
  StoreAccessor &operator=(std::span<const T> buffer)
  /* re-interpret range of memory representing vector of type T as
  a range of memory representing uint8_ts, then write to store */
  {
    return operator=(std::span<const uint8_t>(
        reinterpret_cast<const uint8_t *>(buffer.data()),
        buffer.size() * sizeof(T)));
  }
};

class FSStore
/* A file system (with root in 'basedir' directory) obeying Zarr
version 2 requirements for a Store. Store contins a series
of key, values where values may be data arrays or groups in the store.
data for a given key is written to the store via the functions
in StoreAccessor */
{
private:
  const std::filesystem::path basedir;

public:
  FSStore(std::filesystem::path basedir) : basedir(basedir)
  {
    // initialize a zarr group (i.e. dataset)
    const std::string zarr_format("2"); // storage spec. version 2
    const std::string zgroupjson("{\"zarr_format\": " +
                                   zarr_format + "}");
    (*this)[".zgroup"] = zgroupjson;

    // global metadata (optional)
    (*this)[".zattrs"] = "{\"creator\": \"Clara Bayley\", "
                         "\"title\": \"store for output of coupled SDM\"}";
  }

  StoreAccessor<FSStore> operator[](std::string_view key)
  {
    return {*this, key};
  }

  bool write(std::string_view key, std::span<const uint8_t> buffer);
  /* write function called by StoreAccessor once data has been
  converted into a vector of unsigned integer types */
};

#endif // FSSTORE_HPP