/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: store_accessor.hpp
 * Project: zarr
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors: Tobias Kölling (TB)
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Template class for converting types into vectors of single bytes to
 * write to some kind of memory store under a given key.
 */

#ifndef LIBS_ZARR_STORE_ACCESSOR_HPP_
#define LIBS_ZARR_STORE_ACCESSOR_HPP_

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <span>
#include <string_view>

#include "../kokkosaliases.hpp"

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
  Store& store;         /**< Reference to the store object. */
  std::string_view key; /**< The key under which data will be stored in the store. */

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
  StoreAccessor& operator=(const std::span<const T> buffer) {
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
  StoreAccessor& operator=(const Kokkos::View<T*, HostSpace> buffer) {
    return operator=(std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(buffer.data()),
                                              buffer.extent(0) * sizeof(T)));
  }
};

#endif  // LIBS_ZARR_STORE_ACCESSOR_HPP_
