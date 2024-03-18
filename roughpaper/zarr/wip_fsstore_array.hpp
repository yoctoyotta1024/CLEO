/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: wip_fsstore_array.hpp
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

#ifndef ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_
#define ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_DualView.hpp>

#include "./fsstore.hpp"

using HostSpace = Kokkos::DefaultHostExecutionSpace;     // TODO(CB) (re-)move definitions
using kkpair_size_t = Kokkos::pair<size_t, size_t>;      // TODO(CB) (re-)move definitions

/* converts vector of strings, e.g. for names of dimensions, into a single list
written as a string */
inline std::string vecstr_to_string(const std::vector<std::string> &dims) {
  auto dims_str = std::string{ "[" };
  for (const auto& d : dims) { dims_str += "\"" + d + "\","; }
  dims_str.pop_back();    // delete last ","
  dims_str += "]";
  return dims_str;
}

/* converts vector of integers, e.g. for shape of chunks and array in zarr_metadata, into a single
list written as a string */
inline std::string vec_to_string(const std::vector<size_t> &vals) {
  auto vals_str = std::string{ "[" };
  for (const auto& v : vals) { vals_str += std::to_string(v) + ", "; }
  vals_str.erase(vals_str.size() - 2);    // delete last ", "
  vals_str += "]";
  return vals_str;
}

struct ChunkWriter {
 private:
  std::vector<size_t> chunkshape;
  /**< shape of chunks along each dimension (constant after construction) */
  std::vector<size_t> reducedarray_nchunks;
  /**< number chunks of array along all but outermost dimension (constant after construction) */
  std::vector<size_t> arrayshape;
  /**< number of elements in array along each dimension in store */
  size_t nchunks;
  /**< total number of chunks written in store */

  /* converts vector of integers for label of chunk along each dimension of array
  into a string to use to name a chunk in the store */
  std::string chunk_label() {
    auto chunk_num = std::vector<size_t>(chunkshape.size(), 0);
    chunk_num.at(0) = nchunks / vec_product(reducedarray_nchunks);

    for (size_t aa = 1; aa < chunkshape.size(); ++aa) {
      chunk_num.at(aa) = (nchunks/ vec_product(reducedarray_nchunks, aa)) %
        reducedarray_nchunks.at(aa - 1);
    }

    auto chunk_lab = std::string{ "" };
    for (const auto& c : chunk_num) { chunk_lab += std::to_string(c) + "."; }
    chunk_lab.pop_back();   // delete last "."

    return chunk_lab;
  }

  /* increment shape of outermost dimension of N-Dimensional array and update the array's metadata
  .zarray json correspondingly. Function only updates shape and metadata when the number of chunks
  (nchunks) indicates the reduced array shape is complete, i.e. the array along all but its
  outermost dimension is full of data elements */
  void update_arrayshape(FSStore& store, const std::string_view name,
    const std::string_view partial_metadata, const size_t shape_increment) {
    if (nchunks % vec_product(reducedarray_nchunks) == 0) {
      arrayshape.at(0) += shape_increment;   // increase in shape of outermost dimension
      write_zarray_json(store, name, zarr_metadata(partial_metadata));   // update metadata
    }
  }

 public:
  ChunkWriter(const std::vector<size_t>& chunkshape, const std::vector<size_t>& reduced_arrayshape)
    : chunkshape(chunkshape), arrayshape(chunkshape.size(), 0), nchunks(0) {
    /* number of dimensions for number of chunks must match number of dimensions of array */
    assert((chunkshape.size() == arrayshape.size()));

    /* number of dimensions of reduced array is 1 less than actual array */
    assert((reduced_arrayshape.size() + 1 == chunkshape.size()) &&
      "reduced array 1 less dimension than array (excludes outermost (0th) dimension");

    /* set shape of array and number of chunks along all but array's outermost dimension given
    the shape of each chunk and expected shape of final array along those dimensions */
    for (size_t aa = 1; aa < chunkshape.size(); ++aa) {
      /* Along all but outermost (0th) dimension, the length of a chunk must be completely
      divisible by the array's expected final length along that dimension in order to ensure
      good chunking */
      assert((reduced_arrayshape.at(aa - 1) % chunkshape.at(aa) == 0) &&
        "along all but outermost dimension, arrayshape must be completely divisible by chunkshape");
      /* reducedarray_nchunks = number of chunks along all but outermost dimension of array */
      reducedarray_nchunks.push_back(reduced_arrayshape.at(aa - 1) / chunkshape.at(aa));

      /* set array shape along all but outermost dimension to number of elements given by
      the number and shape of chunks along that dimension */
      arrayshape.at(aa) = reduced_arrayshape.at(aa - 1);
    }
  }

  std::vector<size_t> get_chunkshape() {
    return chunkshape;
  }

  size_t get_reduced_arraysize() {
    return vec_product(arrayshape, 1);
  }

  /* make string of metadata for array in zarr store */
  std::string zarr_metadata(const std::string_view partial_metadata) {
    const auto metadata = std::string(
      "{\n"
      "  \"shape\": " +
      vec_to_string(arrayshape) +
      ",\n" +
      std::string(partial_metadata) +
      "\n}");

    return metadata;
  }

  template <typename T>
  void write_chunk(FSStore& store, const std::string_view name,
    const std::string_view partial_metadata, Buffer<T>& buffer, const size_t shape_increment) {
    buffer.write_buffer_to_chunk(store, name, chunk_label());
    update_arrayshape(store, name, partial_metadata, shape_increment);
    ++nchunks;
  }

  template <typename T>
  void write_chunk(FSStore& store, const std::string_view name,
    const std::string_view partial_metadata, const Buffer<T>::subviewh_buffer h_data_chunk,
    const size_t shape_increment) {
    store[std::string(name) + '/' + chunk_label()].operator=<T>(h_data_chunk);
    update_arrayshape(store, name, partial_metadata, shape_increment);
    ++nchunks;
  }
};

template <typename T>
class FSStoreArrayViaBuffer {
 private:
  using viewh_buffer = Buffer<T>::viewh_buffer;
  using subviewh_buffer = Buffer<T>::subviewh_buffer;
  FSStore& store;                  // file system store satisfying zarr store specificaiton v2
  ChunkWriter chunks;              // information about chunks written in FSStore array
  Buffer<T> buffer;                // buffer for holding data before writing chunks to FSStore array
  std::string_view name;           // name to call variable being stored
  std::string partial_metadata;    // metadata excluding shape required for zarr array

  subviewh_buffer write_chunks_to_store(const subviewh_buffer h_data) {
    // write buffer to chunk if it's full
    if (buffer.get_space() == 0) {
      chunks.write_chunk<T>(store, name, partial_metadata, buffer, chunks.get_chunkshape().at(0));
    }

    // write whole chunks of h_data_remaining
    const auto nchunks_data = size_t{ h_data.extent(0) / buffer.get_chunksize() };
    for (size_t nn = 0; nn < nchunks_data; ++nn) {
      const auto csz = buffer.get_chunksize();
      const auto refs = kkpair_size_t({ nn * csz, (nn + 1) * csz });
      chunks.write_chunk<T>(store, name, partial_metadata, Kokkos::subview(h_data, refs),
        chunks.get_chunkshape().at(0));
    }

    // return remainder of data not written to chunks
    const auto n_to_chunks = nchunks_data * buffer.get_chunksize();
    const auto refs = kkpair_size_t({ n_to_chunks, h_data.extent(0) });
    return Kokkos::subview(h_data, refs);
  }

 public:
  /**
  * @brief Writes a zarr array to a specified file system storafe via a buffer.
  *
  * Initializes an empty array in the provided FSStore in order to writes chunks of array to the
  * store via a buffer. The assertions in this constructor ensure chunks are an appropriate size and
  * shape for the array such that the final array dimensions are exactly integer multiples of its
  * chunks along all but outermost (0th) dimension.
  *
  * @param store The FSStore where the array will be stored.
  * @param chunkshape The shape of individual data chunks along each dimension.
  * @param reduced_arrayshape The shape of the array along all but the outermost (0th) dimension.
  * @param name The name of the array.
  * @param units The units of the array's coordinates.
  * @param scale_factor The scale factor of the data.
  * @param dtype The data type stored in the arrays (e.g., "<f8").
  * @param dims The names of each dimension of the array.
  */
  FSStoreArrayViaBuffer(FSStore& store, const std::vector<size_t>& chunkshape,
    const std::string_view name, const std::string_view units, const double scale_factor,
    const std::string_view dtype, const std::vector<std::string>& dims,
    const std::vector<size_t>& reduced_arrayshape = std::vector<size_t>({}))
    : store(store), chunks(chunkshape, reduced_arrayshape), buffer(chunks.get_chunkshape()),
    name(name) {
    /* number of names of dimensions must match number of dimensions of chunks */
    assert((chunkshape.size() == dims.size()) &&
      "number of named dimensions of array must match number dimensinos of chunks");

    /* chunksize according to buffer must match total size of a (shaped) chunk */
    assert((buffer.get_chunksize() == vec_product(chunks.get_chunkshape())) &&
      "buffer's chunksize must be consistent with chunk shape");

    /* make string of zarray metadata for array in zarr store (incomplete because missing shape) */
    const auto order = 'C';        // layout of bytes in each chunk of array in storage ('C' or 'F')
    const auto compressor = std::string{ "null" };   // compression of data when writing to store
    const auto fill_value = std::string{ "null" };   // fill value for empty datapoints in array
    const auto filters = std::string{ "null" };      // codec configurations for compression
    const auto zarr_format = '2';                    // storage spec. version 2

    partial_metadata = std::string(
      "  \"chunks\": " +
      vec_to_string(chunkshape) +
      ",\n"
      "  \"dtype\": \"" +
      std::string(dtype) +  // dtype = datatype stored in arrays e.g. "<f8"
      "\",\n"
      "  \"order\": \"" +
      order +
      "\",\n"
      "  \"compressor\": " +
      compressor +
      ",\n"
      "  \"fill_value\": " +
      fill_value +
      ",\n"
      "  \"filters\": " +
      filters +
      ",\n"
      "  \"zarr_format\": " +
      zarr_format);

    /* make string of zattrs attribute information for array in zarr store */
    const auto arrayattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dims) +                // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(units) + "\"" +    // units of coordinate being stored
      ",\n"
      "  \"scale_factor\": " +
      std::to_string(scale_factor) +        // scale_factor of data
      "\n}");

    write_zattrs_json(store, name, arrayattrs);
    write_zarray_json(store, name, chunks.zarr_metadata(partial_metadata));
  };

  ~FSStoreArrayViaBuffer() {
    /* write buffer to chunk if it isn't empty */
    if (buffer.get_fill() > 0) {
      const auto reduced_arraysize = chunks.get_reduced_arraysize();   // excluding outer dimension
      if (buffer.get_fill() % reduced_arraysize != 0) {
        const auto warning = std::string_view("WARNING: number of data elements in the buffer"
        " should be completely divisible by the number of elements in the array excluding its"
        " outermost dimension.\n         Some data in this array may be ignored or filled with"
        " null / nan fill value.\n");
        std::cout << warning;
      }
      const auto shape_increment = buffer.get_fill() / reduced_arraysize;
      chunks.write_chunk<T>(store, name, partial_metadata, buffer, shape_increment);
    }
  };

  void write_data_to_zarr_array(const viewh_buffer h_data) {
    auto h_data_rem = buffer.copy_to_buffer(h_data);

    h_data_rem = write_chunks_to_store(h_data_rem);

    h_data_rem = buffer.copy_to_buffer(h_data_rem);

    assert((h_data_rem.extent(0) == 0) && "there is leftover data remaining after writing array");
  };
};

#endif    // ROUGHPAPER_ZARR_WIP_FSSTORE_ARRAY_HPP_
