/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: chunks.hpp
 * Project: zarr
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 20th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Class to manage and write chunks of data to an array in a given memory store.
 */

#ifndef ROUGHPAPER_ZARR_CHUNKS_HPP_
#define ROUGHPAPER_ZARR_CHUNKS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <cassert>
#include <string>
#include <string_view>
#include <vector>

#include "./buffer.hpp"

/**
 * @brief Calculates the product of all elements in a vector of size_t numbers.
 *
 * @param vec The vector of size_t numbers.
 * @return The product of all the elements in the vector.
 */
inline size_t vec_product(const std::vector<size_t>& vec) {
  auto value = size_t{1};
  for (const auto& v : vec) {
    value *= v;
  }
  return value;
}

/**
 * @brief Calculates the product of elements in a vector of size_t numbers
 * starting from the aa'th index of the vector.
 *
 * @param vec The vector of size_t numbers.
 * @param aa The starting index from which to calculate the product.
 * @return The product of elements from the aa'th index in the vector.
 */
inline size_t vec_product(const std::vector<size_t>& vec, const size_t aa) {
  auto value = size_t{1};
  for (auto it = vec.begin() + aa; it != vec.end(); ++it) {
    value *= *it;
  }
  return value;
}

/**
 * @brief A class template for managing and writing chunks of an array.
 *
 * This class provides functionality for writing chunks of an array to a store.
 *
 */
class Chunks {
 private:
  std::vector<size_t> chunkshape; /**< Shape of chunks along each dimension (constant) */
  std::vector<size_t> reducedarray_nchunks;
  /**< Number chunks of array along all but outermost dimension of array (constant) */

  /**
   * @brief Create label for a chunk given current number of chunks written to array.
   *
   * This function creates a vector of integers for the number of a chunk along each dimension of
   * an array given the chunk is the n'th chunk to be written to the store (starting at n=0). The
   * vector is then converted into a string which can be used to label the chunk.
   *
   * @param chunk_num The number of the chunk to write to the array.
   * @return A string representing the label of the current chunk to write.
   */
  std::string chunk_label(const size_t chunk_num) const {
    auto chunk_labnums = std::vector<size_t>(chunkshape.size(), 0);
    chunk_labnums.at(0) = chunk_num / vec_product(reducedarray_nchunks);

    for (size_t aa = 1; aa < chunkshape.size(); ++aa) {
      chunk_labnums.at(aa) =
          (chunk_num / vec_product(reducedarray_nchunks, aa)) % reducedarray_nchunks.at(aa - 1);
    }

    auto chunk_lab = std::string{""};
    for (const auto& c : chunk_labnums) {
      chunk_lab += std::to_string(c) + ".";
    }
    chunk_lab.pop_back();  // delete last "."

    return chunk_lab;
  }

 public:
  /**
   * @brief Constructor for the Chunks class.
   *
   * Initializes the Chunks with the provided chunk shape and reduced array shape. Reduced
   * array shape is the shape of the array along all but the outermost dimensions of the array.
   *
   * @param chunkshape The shape of chunks along each dimension.
   * @param reduced_arrayshape The shape of the reduced array along each dimension.
   */
  Chunks(const std::vector<size_t>& chunkshape, const std::vector<size_t>& reduced_arrayshape)
      : chunkshape(chunkshape), reducedarray_nchunks(chunkshape.size() - 1, 0) {
    /* number of dimensions of reduced array is 1 less than actual array ( = array's chunks) */
    assert((reduced_arrayshape.size() == chunkshape.size() - 1) &&
           "reduced array 1 less dimension than array (excludes outermost (0th) dimension");

    /* set number of chunks along all but array's outermost dimension given
    the shape of each chunk and expected shape of final array along those dimensions */
    for (size_t aa = 1; aa < chunkshape.size(); ++aa) {
      /* Assert the chunk size is completely divisible by the array's expected size along that
      dimension in order to ensure good chunking */
      assert((reduced_arrayshape.at(aa - 1) % chunkshape.at(aa) == 0) &&
             "along all but outermost dimension, arrayshape must be completely divisible by "
             "chunkshape");
      /* reducedarray_nchunks = number of chunks along all but outermost dimension of array */
      reducedarray_nchunks.push_back(reduced_arrayshape.at(aa - 1) / chunkshape.at(aa));
    }
  }

  /**
   * @brief Gets the shape of a chunk.
   *
   * @return A vector containing the shape (number of data elements) of a chunk
   * along each dimension.
   */
  std::vector<size_t> get_chunkshape() const { return chunkshape; }

  /**
   * @brief Gets the number of chunks of the reduced array.
   *
   * @return A vector containing the number of chunks of an array along its dimensions except for
   * its outermost one.
   */
  std::vector<size_t> get_reducedarray_nchunks() const { return reducedarray_nchunks; }

  /**
   * @brief Writes a chunk to the store and increments the total number of chunks written.
   *
   * This function writes the data held in a buffer in the specified store to a chunk identified by
   * "chunk_label" of an array called "name" given the number of chunks of the array already
   * existing. After writing the chunk, the total number of chunks is incremented.
   *
   * @tparam Store The type of the store.
   * @tparam T The type of the data elements stored in the buffer.
   * @param store Reference to the store where the chunk will be written.
   * @param name Name of the array in the store where the chunk will be written.
   * @param totnchunks The total number of chunks of the array already written.
   * @param buffer The buffer containing the data to be written to the chunk.
   * @return The updated total number of chunks after writing.
   */
  template <typename Store, typename T>
  size_t write_chunk(Store& store, const std::string_view name, const size_t totnchunks,
                     Buffer<T>& buffer) const {
    buffer.write_buffer_to_chunk(store, name, chunk_label(totnchunks));
    return totnchunks + 1;
  }

  /**
   * @brief Writes a chunk to the store and increments the total number of chunks written.
   *
   * This function writes the data stored in the Kokkos view (in host memory) in the specified store
   * to a chunk identified by "chunk_label" of an array called "name" given the number of chunks of
   * the array already existing. After writing the chunk, the total number of chunks is incremented.
   *
   * @tparam Store The type of the store.
   * @tparam T The type of the data elements stored in the buffer.
   * @param store Reference to the store where the chunk will be written.
   * @param name Name of the array in the store where the chunk will be written.
   * @param totnchunks The total number of chunks of the array already written.
   * @param h_data_chunk The view containing the data in host memory to be written to the chunk.
   * @return The updated total number of chunks after writing.
   */
  template <typename Store, typename T>
  size_t write_chunk(Store& store, const std::string_view name, const size_t totnchunks,
                     const Buffer<T>::subviewh_buffer h_data_chunk) const {
    store[std::string(name) + '/' + chunk_label(totnchunks)].operator= <T>(h_data_chunk);
    return totnchunks + 1;
  }
};

#endif  // ROUGHPAPER_ZARR_CHUNKS_HPP_