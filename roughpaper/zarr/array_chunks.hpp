/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: array_chunks.hpp
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
 * Class to manage and write chunks of data to an array in a given memory store.
 */


#ifndef ROUGHPAPER_ZARR_ARRAY_CHUNKS_HPP_
#define ROUGHPAPER_ZARR_ARRAY_CHUNKS_HPP_

#include <vector>
#include <string>
#include <string_view>

#include "./buffer.hpp"

/**
 * @brief A class template for managing and writing chunks of an array.
 *
 * This class provides functionality for writing chunks of an array to a store.
 *
 * @tparam T The type of data elements stored in the buffer.
 */
template <typename T>
class ArrayChunks {
 private:
  std::vector<size_t> chunkshape;   /**< Shape of chunks along each dimension */
  std::vector<size_t> reducedarray_nchunks;
  /**< Number chunks of array along all but outermost dimension of array */

  /* converts vector of integers for label of chunk along each dimension of array into a string
  to use to name a chunk in the store */
  /**
   * @brief Create label for a chunk given current number of chunks written to array.
   *
   * This function creates and converts a vector of integers representing the label of a
   * chunk along each dimension of an array into a string which can be used to name the current
   * chunk that is next to be written to the store.
   *
   * @return A string representing the label of the current chunk to write.
   */
  std::string chunk_label(const size_t totnchunks) {
    auto chunk_num = std::vector<size_t>(chunkshape.size(), 0);
    chunk_num.at(0) = totnchunks / vec_product(reducedarray_nchunks);

    for (size_t aa = 1; aa < chunkshape.size(); ++aa) {
      chunk_num.at(aa) = (totnchunks/ vec_product(reducedarray_nchunks, aa)) %
        reducedarray_nchunks.at(aa - 1);
    }

    auto chunk_lab = std::string{ "" };
    for (const auto& c : chunk_num) { chunk_lab += std::to_string(c) + "."; }
    chunk_lab.pop_back();   // delete last "."

    return chunk_lab;
  }

 public:
  /**
   * @brief Constructor for the ArrayChunks class.
   *
   * Initializes the ArrayChunks with the provided chunk shape and reduced array shape. Reduced
   * array shape is the shape of the array along all but the outermost dimensions of the array.
   *
   * @param chunkshape The shape of chunks along each dimension.
   * @param reduced_arrayshape The shape of the reduced array along each dimension.
   */
  ChunkWriter(const std::vector<size_t>& chunkshape, const std::vector<size_t>& reduced_arrayshape)
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
        "along all but outermost dimension, arrayshape must be completely divisible by chunkshape");
      /* reducedarray_nchunks = number of chunks along all but outermost dimension of array */
      reducedarray_nchunks.push_back(reduced_arrayshape.at(aa - 1) / chunkshape.at(aa));
    }
  }

  std::vector<size_t> get_chunkshape() {
    return chunkshape;
  }

  template <typename Store, typename T>
  size_t write_chunk(Store& store, const std::string_view name, const size_t totnchunks,
    Buffer<T>& buffer) {
    buffer.write_buffer_to_chunk(store, name, chunk_label(totnchunks));
    return ++totnchunks;
  }

  template <typename Store, typename T>
  size_t write_chunk(Store& store, const std::string_view name, const size_t totnchunks,
    const Buffer<T>::subviewh_buffer h_data_chunk) {
    store[std::string(name) + '/' + chunk_label(totnchunks)].operator=<T>(h_data_chunk);
    return ++totnchunks;
  }
};

#endif   // ROUGHPAPER_ZARR_ARRAY_CHUNKS_HPP_
