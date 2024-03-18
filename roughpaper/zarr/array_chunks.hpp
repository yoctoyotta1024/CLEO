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
 * Class to count and write chunks of data to an array in a given memory store.
 */


#ifndef ROUGHPAPER_ZARR_ARRAY_CHUNKS_HPP_
#define ROUGHPAPER_ZARR_ARRAY_CHUNKS_HPP_

#include <vector>
#include <string>
#include <string_view>

#include "./buffer.hpp"

template <typename T>
class ArrayChunks {
 private:
  size_t totnchunks;
  /**< total number of chunks written in store */
  std::vector<size_t> chunkshape;
  /**< shape of chunks along each dimension (constant) */
  std::vector<size_t> reducedarray_nchunks;
  /**< number chunks of array along all but outermost dimension of array (constant) */

  /* converts vector of integers for label of chunk along each dimension of array into a string
  to use to name a chunk in the store */
  std::string chunk_label() {
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
  ChunkWriter(const std::vector<size_t>& chunkshape, const std::vector<size_t>& reduced_arrayshape)
    : totnchunks(0), chunkshape(chunkshape), reducedarray_nchunks(chunkshape.size() - 1, 0) {

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
  void write_chunk(Store& store, const std::string_view name, Buffer<T>& buffer) {
    buffer.write_buffer_to_chunk(store, name, chunk_label());
    ++totnchunks;
  }

  template <typename Store, typename T>
  void write_chunk(Store& store, const std::string_view name,
    const Buffer<T>::subviewh_buffer h_data_chunk) {
    store[std::string(name) + '/' + chunk_label()].operator=<T>(h_data_chunk);
    ++totnchunks;
  }
};

#endif   // ROUGHPAPER_ZARR_ARRAY_CHUNKS_HPP_
