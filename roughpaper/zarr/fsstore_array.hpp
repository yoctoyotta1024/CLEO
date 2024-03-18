/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: fsstore_array.hpp
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

#ifndef ROUGHPAPER_ZARR_FSSTORE_ARRAY_HPP_
#define ROUGHPAPER_ZARR_FSSTORE_ARRAY_HPP_

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

/* returns product of a vector of size_t numbers */
inline size_t vec_product(const std::vector<size_t>& vec) {
  auto value = size_t{1};
  for (const auto& v : vec) {
    value *= v;
  }
  return value;
}

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

template <typename T>
struct Buffer {
 public:
  using viewh_buffer = Kokkos::View<T*, HostSpace::memory_space>;     // view of buffer type on host
  using subviewh_buffer = Kokkos::Subview<viewh_buffer, kkpair_size_t>;  // subview of host view

 private:
  size_t chunksize;                        // total chunk size = product of shape of chunks
  size_t fill;                             // number of elements of buffer currently filled
  viewh_buffer buffer;                     // view for buffer in host memory

  /* parallel loop on host to fill buffer with nan (numerical limit) values */
  void reset_buffer() {
    Kokkos::parallel_for(
      "init_buffer", Kokkos::RangePolicy<HostSpace>(0, chunksize),
      KOKKOS_CLASS_LAMBDA(const size_t & jj) {
      buffer(jj) = std::numeric_limits<T>::max();
    });
    fill = 0;
  }

  /* parallel loop on host to fill buffer from start of empty spaces (i.e. from index "fill")
  with "n_to_copy" elements from data */
  void copy_ndata_to_buffer(const size_t n_to_copy, const viewh_buffer h_data) {
    Kokkos::parallel_for(
      "copy_ndata_to_buffer", Kokkos::RangePolicy<HostSpace>(fill, fill + n_to_copy),
      KOKKOS_CLASS_LAMBDA(const size_t & jj) {
      buffer(jj) = h_data(jj);
    });
    fill += n_to_copy;
  }

 public:
  explicit Buffer(const std::vector<size_t>& chunkshape) : chunksize(vec_product(chunkshape)),
    fill(0), buffer("buffer", chunksize) {
    reset_buffer();
  }

  size_t get_chunksize() {
    return chunksize;
  }

  size_t get_fill() {
    return fill;
  }

  /* returns number of spaces in buffer currently not filled with data */
  size_t get_space() {
    return chunksize - fill;
  }

  /* copies as many as possible elements of data to buffer until either all the data is written to
  the buffer, or all the spaces in the buffer are filled. Returns view of remaining data not copied
  to buffer (empty if all the data is copied). */
  subviewh_buffer copy_to_buffer(const viewh_buffer h_data) {
    // number of elements of data to copy to buffer
    const auto n_to_copy = size_t{ std::min(get_space(), h_data.extent(0)) };

    // copy "n_to_copy" number of elements of data to buffer
    copy_ndata_to_buffer(n_to_copy, h_data);

    // return remainder of data not copied to buffer
    const auto refs = kkpair_size_t({ n_to_copy, h_data.extent(0) });
    return Kokkos::subview(h_data, refs);
  }

  /* write out data from buffer to chunk called "chunk_str" in an array called "name" in a (zarr)
  file system store. Then reset buffer. */
  void write_buffer_to_chunk(FSStore& store, std::string_view name, const std::string &chunk_str) {
    store[std::string(name) + '/' + chunk_str].operator=<T>(buffer);
    reset_buffer();
  }
};

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

  /* converts vector of integers for chunkcount into string to use to name a chunk */
  std::string chunks_string() {
    std::vector<size_t> chunkcount;   // TODO(CB) # of chunks along each dimension written in store
    auto chunk_str = std::string{ "" };
    for (const auto& c : chunkcount) { chunk_str += std::to_string(c) + "."; }
    chunk_str.pop_back();   // delete last "."

    return chunk_str;
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
    assert((reduced_arrayshape.size() + 1 == arrayshape.size()) &&
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
    const std::string_view partial_metadata, Buffer<T>& buffer, const std::vector<size_t>& shape) {
    buffer.write_buffer_to_chunk(store, name, chunks_string());
    update_arrayshape(store, name, partial_metadata, shape);
    ++nchunks;
  }

  template <typename T>
  void write_chunk(FSStore& store, const std::string_view name,
    const std::string_view partial_metadata, const Buffer<T>::subviewh_buffer h_data_chunk,
    const std::vector<size_t>& shape) {
    const auto chunk_str = chunks_string();
    store[std::string(name) + '/' + chunk_str].operator=<T>(h_data_chunk);
    update_arrayshape(store, name, partial_metadata, shape);
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
      chunks.write_chunk<T>(store, name, partial_metadata, buffer, chunks.get_chunkshape());
    }

    // write whole chunks of h_data_remaining
    const auto nchunks_data = size_t{ h_data.extent(0) / buffer.get_chunksize() };
    for (size_t nn = 0; nn < nchunks_data; ++nn) {
      const auto csz = buffer.get_chunksize();
      const auto refs = kkpair_size_t({ nn * csz, (nn + 1) * csz });
      chunks.write_chunk<T>(store, name, partial_metadata, Kokkos::subview(h_data, refs),
        chunks.get_chunkshape());
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
    // write buffer to chunk if it isn't empty
    if (buffer.get_fill() > 0) {
      const auto chunkshape = chunks.get_chunkshape();
      auto shape = std::vector<size_t>(chunkshape.begin() + 1, chunkshape.end());
      assert((buffer.get_fill() % vec_product(shape) == 0) &&
        "data in buffer should be completely divisible by reduced chunkshape");
      shape.insert(shape.begin(), buffer.get_fill() / vec_product(shape));
      chunks.write_chunk<T>(store, name, partial_metadata, buffer, shape);
    }
  };

  void write_data_to_zarr_array(const viewh_buffer h_data) {
    auto h_data_rem = buffer.copy_to_buffer(h_data);

    h_data_rem = write_chunks_to_store(h_data_rem);

    h_data_rem = buffer.copy_to_buffer(h_data_rem);

    assert((h_data_rem.extent(0) == 0) && "there is leftover data remaining after writing array");
  };
};

#endif    // ROUGHPAPER_ZARR_FSSTORE_ARRAY_HPP_
