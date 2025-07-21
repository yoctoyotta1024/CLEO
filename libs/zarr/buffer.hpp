/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: buffer.hpp
 * Project: zarr
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Class for a buffer used by a ZarrArray to acculuate data and then write it into a store
 */

#ifndef LIBS_ZARR_BUFFER_HPP_
#define LIBS_ZARR_BUFFER_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <algorithm>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "../kokkosaliases.hpp"

/**
 * @brief A class template for managing a buffer of elements of data type T.
 *
 * This class provides functionality for initializing a buffer, copying elements of data into it
 * and writing the buffer to a store.
 *
 * @tparam The type of the store object used by the buffer.
 * @tparam T The type of elements stored in the buffer.
 */
template <typename T>
struct Buffer {
 public:
  using viewh_buffer = Kokkos::View<T*, HostSpace>; /**< View of buffer type on host */
  using subviewh_buffer = Kokkos::Subview<viewh_buffer, kkpair_size_t>; /**< Subview of host view */
  using mirrorviewd_buffer = Kokkos::View<T*, HostSpace::array_layout,
                                          ExecSpace>; /**< mirror view of buffer view on device */

 private:
  const size_t chunksize; /**< Total chunk size = product of shape of chunks */
  size_t fill;            /**< Number of elements of buffer currently filled */
  viewh_buffer buffer;    /**< View for buffer in host memory */

  /**
   * @brief Parallel loop on host to fill buffer with NaN (numerical limit).
   */
  void reset_buffer() {
    auto buffer_ = buffer;  // Copy of view for lambda functions using buffer

    Kokkos::parallel_for(
        "reset_buffer", Kokkos::RangePolicy<HostSpace>(0, chunksize),
        [buffer_](const size_t& jj) { buffer_(jj) = std::numeric_limits<T>::max(); });
    fill = 0;
  }

  /**
   * @brief Parallel loop on host to fill buffer with data elements
   *
   * Parallel loop on host to fill buffer from index "fill" (i.e. start of empty spaces)
   * with "n_to_copy" elements from view of data.
   *
   * @param n_to_copy maximum number of elements to copy to the buffer.
   * @param h_data View containing the data to copy.
   */
  void copy_ndata_to_buffer(const size_t n_to_copy, const viewh_buffer h_data) {
    const auto refs_d = kkpair_size_t({0, n_to_copy});
    const auto source = Kokkos::subview(h_data, refs_d);  // data to copy into buffer

    const auto refs_b = kkpair_size_t({fill, fill + n_to_copy});
    const auto dest = Kokkos::subview(buffer, refs_b);  // space in buffer to paste into

    Kokkos::deep_copy(dest, source);

    fill += n_to_copy;
  }

 public:
  /**
   * @brief Constructor for the Buffer class.
   *
   * Initializes the buffer with size of given chunkshape.
   *
   * @param chunksize number of elements of data in 1 chunk of an array.
   */
  explicit Buffer(const size_t chunksize)
      : chunksize(chunksize), fill(0), buffer("buffer", chunksize) {
    reset_buffer();
  }

  /**
   * @brief Gets the total chunk size of the buffer.
   *
   * @return The total chunk size.
   */
  size_t get_chunksize() const { return chunksize; }

  /**
   * @brief Gets the number of elements currently in the buffer.
   *
   * @return The number of elements of buffer filled.
   */
  size_t get_fill() const { return fill; }

  /**
   * @brief Returns the number of empty spaces in the buffer.
   *
   * @return The number of spaces in the buffer currently not filled with data.
   */
  size_t get_space() const { return chunksize - fill; }

  /**
   * @brief Copies as many elements as possible from data to buffer.
   *
   * Copies "n_to_copy" elements from view of data such that either all the data is copied
   * to the buffer or all the spaces in the buffer are filled. Returns a view of remaining
   * data not copied to the buffer which is empty if all the data has been copied.
   *
   * @param h_data View containing the data to copy.
   * @return Subview containing the remaining data not copied to the buffer.
   */
  subviewh_buffer copy_to_buffer(const viewh_buffer h_data) {
    const auto n_to_copy = size_t{std::min(get_space(), h_data.extent(0))};

    copy_ndata_to_buffer(n_to_copy, h_data);

    const auto refs = kkpair_size_t({n_to_copy, h_data.extent(0)});  // indexes of remaining data
    return Kokkos::subview(h_data, refs);
  }

  /**
   * @brief Copies maximum of 1 element of data to buffer.
   *
   * Assert that there is space in the buffer, then copy 1 element of data to
   * the buffer.
   *
   * @param data Data element to copy.
   */
  void copy_to_buffer(const T data) {
    assert((get_space() > 0) && "buffer must have space to copy element");
    buffer(fill) = data;
    ++fill;
  }

  /**
   * @brief Writes data from buffer to a chunk in a store.
   *
   * Writes data from buffer to a chunk specified by "chunk_label" of an array
   * called "name" in a memory store. Then resets the buffer.
   *
   * @tparam Store The type of the memory store.
   * @param store Reference to the store object.
   * @param name Name of the array in the store.
   * @param chunk_label Name of the chunk of the array to write in the store.
   */
  template <typename Store>
  void write_buffer_to_chunk(Store& store, std::string_view name, const std::string& chunk_label) {
    store[std::string(name) + '/' + chunk_label] = buffer;
    reset_buffer();
  }
};

#endif  // LIBS_ZARR_BUFFER_HPP_
