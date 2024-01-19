/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: fsstore.hpp
 * Project: zarr
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 23rd October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * objects that can be used as stores obeying the
 * zarr storage specification version 2 (e.g. see FSStore)
 * https://zarr.readthedocs.io/en/stable/spec/v2.html
 */

#ifndef LIBS_ZARR_FSSTORE_HPP_
#define LIBS_ZARR_FSSTORE_HPP_

#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <string>
#include <string_view>

/* functions for converting types (e.g. vectors of
unsigned integers or doubles) into vectors of single bytes to
write to store under a given key. Store can be anything that
satisfies the zarr storage specifcaiton version 2 */
template <typename Store>
struct StoreAccessor {
  Store &store;
  std::string_view key;

  /* write range of memory representing uint8_ts to store */
  StoreAccessor &operator=(std::span<const uint8_t> buffer) {
    store.write(key, buffer);
    return *this;
  }

  /* reinterpret range of memory representing string as
  a range of memory representing uint8_ts, then write to store */
  StoreAccessor &operator=(std::string_view buffer) {
    return operator=(
        std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(buffer.data()), buffer.size()));
  }

  /* re-interpret range of memory representing vector of type T as
  a range of memory representing uint8_ts, then write to store */
  template <typename T>
  StoreAccessor &operator=(std::span<const T> buffer) {
    return operator=(std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(buffer.data()),
                                              buffer.size() * sizeof(T)));
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
    const std::string zarr_format("2");  // storage spec. version 2
    const std::string zgroupjson("{\"zarr_format\": " + zarr_format + "}");
    (*this)[".zgroup"] = zgroupjson;

    // global metadata (optional)
    (*this)[".zattrs"] =
        "{\"creator\": \"Clara Bayley\", "
        "\"title\": \"store for output of coupled SDM\"}";
  }

  StoreAccessor<FSStore> operator[](const std::string_view key) { return {*this, key}; }

  /* write function called by StoreAccessor once data has been
  converted into a vector of unsigned integer types */
  bool write(const std::string_view key, const std::span<const uint8_t> buffer);
};

#endif  // LIBS_ZARR_FSSTORE_HPP_
