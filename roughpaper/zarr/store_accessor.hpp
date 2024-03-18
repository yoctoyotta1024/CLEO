/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: store_accessor.hpp
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
 * Template class for converting types into vectors of single bytes to
 * write to some kind of memory store under a given key.
 */


#ifndef ROUGHPAPER_ZARR_STORE_ACCESSOR_HPP_
#define ROUGHPAPER_ZARR_STORE_ACCESSOR_HPP_

#include <span>
#include <string_view>

#include <Kokkos_Core.hpp>

/**
 * @brief A template class for converting types into vectors of single bytes to write to a "store"
 * under a given key.
 *
 * This class provides functions for converting various types (e.g., vectors / Kokkos views of
 * unsigned integers or doubles) into spans of single bytes, which can then be written to a "store"
 * under a specified key. The store is any object, but for example it could be a file system
 * satisfying the Zarr storage specification version 2.
 *
 * @tparam Store The type of the store object.
 */
template <typename Store>
struct StoreAccessor {
  using HostSpace = Kokkos::DefaultHostExecutionSpace;    // TODO(CB) (re-)move definitions

  Store& store;   ///< Reference to the store object.
  std::string_view key;   ///< The key under which data will be stored in the store.

  /**
   * @brief Write a range of memory representing unsigned bytes (uint8_t) to the store.
   *
   * @param buffer A span representing the range of memory containing the unsigned bytes.
   * @return A reference to the current StoreAccessor object.
   */
  StoreAccessor& operator=(std::span<const uint8_t> buffer) {
    store.write(key, buffer);
    return *this;
  }

  /**
   * @brief Convert a C++ string type into a range of memory representing unsigned bytes,
   * then writes the bytes to the store.
   *
   * Reinterprets (i.e. converts) the range of memory occupied by a string as a series of unsigned
   * bytes (uint8_t), and then writes this memory to the store.
   *
   * @param buffer The string to be converted and written to the store.
   * @return A reference to the current StoreAccessor object.
   */
  StoreAccessor& operator=(std::string_view buffer) {
    return operator=(
      std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(buffer.data()), buffer.size()));
  }

  /**
   * @brief Convert a vector of type T into a range of memory representing unsigned bytes
   * then writes the bytes to the store.
   *
   * Reinterprets (i.e. converts) the range of memory occupied by a vector of type T as a series of
   * unsigned bytes (uint8_t), and then writes this memory to the store.
   *
   * @tparam T The type of the elements in the vector.
   * @param buffer The vector to be converted and written to the store.
   * @return A reference to the current StoreAccessor object.
   */
  template <typename T>
  StoreAccessor& operator=(std::span<const T> buffer) {
    return operator=(std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(buffer.data()),
      buffer.size() * sizeof(T)));
  }

  /**
   * @brief Convert a Kokkos view of type T in the host memory space into a range of memory
   * representing unsigned bytes, then writes the bytes to the store.
   *
   * Reinterprets (i.e. converts) the range of memory occupied by the elements of a Kokkos view of
   * type T in the host memory space as a series of unsigned bytes (uint8_t), and then writes this
   * memory to the store.
   *
   * @tparam T The type of the elements in the Kokkos view in the host memory space.
   * @param buffer The Kokkos view to be converted and written to the store.
   * @return A reference to the current StoreAccessor object.
   */
  template <typename T>
  StoreAccessor& operator=(const Kokkos::View<T*, HostSpace::memory_space> buffer) {
    return operator=(std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(buffer.data()),
      buffer.extent(0) * sizeof(T)));
  }
};


/**
 * @brief Write metadata string to a store under a .zarray key.
 *
 * write metadata under .zarray key in store for an array called 'name'. The key and metadata
 * could be anything, but for example .zarray could be a json file in a file system store
 * (see FSStore) for the metadata which must exist in order to decode chunks of an array according
 * to zarr storage specification version 2 (https://zarr.readthedocs.io/en/stable/spec/v2.html),
 *
 * @tparam Store The type of the store object where the metadata will be written.
 * @param store The store object where the metadata will be written.
 * @param name The name under which the .zarray key will be stored in the store.
 * @param metadata The metadata to write for the .zarray key.
 */
template <typename Store>
inline void write_zarray_json(Store& store, std::string_view name, std::string_view metadata) {
  store[std::string(name) + "/.zarray"] = metadata;
}

/**
 * @brief Write attributes string to a store under a .zattrs key.
 *
 * Write some data under .zattrs key in store for an array called 'name'. The key and attrs data
 * could be anything, but for example .zattrs could be a json file in a file system store
 * (see FSStore) for the extra metadata which must exist in order to make xarray and netCDF
 * happy when opening a zarr dataset, e.g. by naming the dimensions of the
 * "{\"_ARRAY_DIMENSIONS\": [\"dimension_name\"]}";.
 *
 * @tparam Store The type of the store object where the metadata will be written.
 * @param store The store object where the metadata will be written.
 * @param name The name under which the .zarray key will be stored in the store.
 * @param metadata The metadata to write for the .zarray key.
 */
template <typename Store>
inline void write_zattrs_json(Store& store, std::string_view name, std::string_view attrs) {
  store[std::string(name) + "/.zattrs"] = attrs;
}

#endif   // ROUGHPAPER_ZARR_STORE_ACCESSOR_HPP_
