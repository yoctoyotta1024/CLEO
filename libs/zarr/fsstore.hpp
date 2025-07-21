/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: fsstore.hpp
 * Project: zarr
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors: Tobias Kölling (TB)
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Class for writing memory in a a file system store under a given key.
 */

#ifndef LIBS_ZARR_FSSTORE_HPP_
#define LIBS_ZARR_FSSTORE_HPP_

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <string_view>

#include "zarr/store_accessor.hpp"

/**
 * @brief A file system store e.g. for Zarr arrays or groups.
 *
 * This class represents a file system store for a series of key-value pairs for example for
 * storing Zarr data arrays or groups. Data for a given key is written to the store via helper
 * functions from the StoreAccessor struct.
 *
 */
class FSStore {
 private:
  const std::filesystem::path basedir; /**< The root directory of the file system store. */

 public:
  /**
   * @brief Constructs an FSStore object with the specified base directory.
   *
   * @param basedir The root directory of the file system store.
   */
  explicit FSStore(const std::filesystem::path basedir) : basedir(basedir) {}

  /**
   * @brief Operator to use a StoreAccessor to write values under a given key.
   *
   * Usage: `operator[y] = x;` writes values x under a key called 'y' using the StoreAccessor.
   *
   * @param key The key for which the StoreAccessor is accessed.
   * @return A StoreAccessor object associated with the specified key.
   */
  StoreAccessor<FSStore> operator[](const std::string_view key) { return {*this, key}; }

  /**
   * @brief Write function called by StoreAccessor to write data to file system storage after the
   * data has been converted into a vector of unsigned integer types.
   *
   * This function can be used by a StoreAccessor object to write data represented as a vector of
   * unsigned integer types to the file system store under the specified key.
   *
   * @param key The key under which the data will be stored in the file system store.
   * @param buffer A span representing the range of memory containing the unsigned bytes to be
   * written.
   * @return True if the write operation is successful, false otherwise.
   */
  bool write(const std::string_view key, const std::span<const uint8_t> buffer) const;
};

#endif  // LIBS_ZARR_FSSTORE_HPP_
