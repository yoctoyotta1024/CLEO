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
 * Last Modified: Thursday 14th March 2024
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

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_DualView.hpp>

#include "./fsstore.hpp"

using dualview_type = Kokkos::DualView<double*>;             // dual view of doubles

using kkpair_size_t = Kokkos::pair<size_t, size_t>;
using subview_type = Kokkos::Subview<dualview_type::t_host, kkpair_size_t>;  // subview of host view

using HostSpace = Kokkos::DefaultHostExecutionSpace;
using viewh_buffer = Kokkos::View<double*, HostSpace::memory_space>;   // view for buffer on host

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

struct Buffer {
 private:
  size_t chunksize;                        // total chunk size = product of shape of chunks
  size_t fill;                             // number of elements of buffer currently filled
  viewh_buffer buffer;

  /* parallel loop on host to fill buffer with nan (numerical limit) values */
  void reset_buffer() {
    Kokkos::parallel_for(
      "init_buffer", Kokkos::RangePolicy<HostSpace>(0, chunksize),
      KOKKOS_CLASS_LAMBDA(const size_t & jj) {
      buffer(jj) = std::numeric_limits<double>::max();
    });
    fill = 0;
  }

  /* parallel loop on host to fill buffer from start of empty spaces (i.e. from index "fill")
  with "n_to_copy" elements from data */
  void copy_ndata_to_buffer(const size_t n_to_copy, const dualview_type::t_host h_data) {
    Kokkos::parallel_for(
      "copy_ndata_to_buffer", Kokkos::RangePolicy<HostSpace>(fill, fill + n_to_copy),
      KOKKOS_CLASS_LAMBDA(const size_t & jj) {
      buffer(jj) = h_data(jj);
    });
    fill += n_to_copy;
  }

 public:
  explicit Buffer(const std::vector<size_t>& chunkshape) : chunksize(1), fill(0),
    buffer("buffer", chunksize) {
    for (const auto& c : chunkshape) { chunksize *= c; }
    Kokkos::resize(buffer, chunksize);
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
  subview_type copy_to_buffer(const dualview_type::t_host h_data) {
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
  void write_buffer_to_chunk(FSStore& store, std::string_view name, std::string_view chunk_str) {
    std::cout << "--> writing buffer to chunk: " << chunk_str << "\n";
    // store[name + "/" + chunk_str].operator= <T>(buffer); // TODO(CB) write buffer chunk
    reset_buffer();
  }
};

struct ChunkWriter {
 private:
  std::vector<size_t> chunkcount;          // number of chunks along each dimension written in store
  std::vector<size_t> arrayshape;          // number of elements in array along each dimension

  /* converts vector of integers for chunkcount into string to use to name a chunk */
  std::string chunkcount_to_string() {
    auto chunk_str = std::string{ "" };
    for (const auto& c : chunkcount) { chunk_str += std::to_string(c) + "."; }
    chunk_str.pop_back();   // delete last "."

    return chunk_str;
  }

  /* update numbers of chunks and shape of array along each dimension */
  void update_chunks_and_shape(const std::vector<size_t>& shape_incre,
    const std::vector<size_t>& chunk_incre) {
    for (size_t aa = 0; aa < chunkcount.size(); ++aa) {
      chunkcount.at(aa) += chunk_incre.at(aa);
      arrayshape.at(aa) += shape_incre.at(aa);
    }
  }

 public:
  explicit ChunkWriter(const size_t ndims) : chunkcount(std::vector<size_t>(ndims, 0)),
    arrayshape(std::vector<size_t>(ndims, 0)) {}

  std::string get_arrayshape_str() {
    return vec_to_string(arrayshape);
  }

  void write_chunk(FSStore& store, std::string_view name, Buffer &buffer) {
    const auto chunk_str = chunkcount_to_string();
    const auto chunk_incre = std::vector<size_t>({1});    // TODO(CB) deal with multi-D chunks
    const auto shape_incre = std::vector<size_t>({buffer.get_fill()});   // TODO(CB) dito
    buffer.write_buffer_to_chunk(store, name, chunk_str);
    update_chunks_and_shape(shape_incre, chunk_incre);
  }

  void write_chunk(FSStore& store, std::string_view name, const subview_type h_data_chunk) {
    const auto chunk_str = chunkcount_to_string();
    const auto chunk_incre = std::vector<size_t>({1});    // TODO(CB) deal with multi-D chunks
    const auto shape_incre = std::vector<size_t>({h_data_chunk.extent(0)});   // TODO(CB) dito
    std::cout << "--> writing h_data to chunk: " << chunk_str << "\n";
    // write_data_to_chunk(store, name, chunk_str, h_data_chunk);   // TODO(CB) write subview chunk
    update_chunks_and_shape(shape_incre, chunk_incre);
  }
};

class FSStoreArrayViaBuffer {
 private:
  FSStore& store;                  // file system store satisfying zarr store specificaiton v2
  Buffer buffer;                   // buffer for holding data before writing chunks to FSStore array
  ChunkWriter chunks;              // information about chunks written in FSStore array
  std::string_view name;           // name to call variable being stored
  std::vector<size_t> chunkshape;  // shape of chunks of array along each dimension
  std::string partial_metadata;    // metadata excluding shape required for zarr array

  /* make string of metadata for array in zarr store */
  std::string zarr_metadata() {
    const auto metadata = std::string(
      "{\n"
      "\"shape\": " +
      chunks.get_arrayshape_str() +
      ",\n" +
      partial_metadata +
      ",\n}");

    return metadata;
  }

  subview_type write_chunks_in_store(const subview_type h_data) {
    // write buffer to chunk if it's full
    if (buffer.get_space() == 0) {
      chunks.write_chunk(store, name, buffer);
    }

    // write whole chunks of h_data_remaining
    const auto nchunks_data = size_t{ h_data.extent(0) / buffer.get_chunksize() };
    std::cout << "nchunks from h_data: " << nchunks_data << "\n";
    for (size_t bb = 0; bb < nchunks_data; ++bb) {
      const auto start = size_t{ bb * buffer.get_chunksize() };
      const auto end = size_t{start + buffer.get_chunksize()};
      const auto refs = kkpair_size_t({ start, end });
      chunks.write_chunk(store, name, Kokkos::subview(h_data, refs));
    }

    // update zarry json with new metadata now chunks have been written
    write_zarray_json(store, name, zarr_metadata());

    // return remainder of data not written to chunks
    const auto n_to_chunks = nchunks_data * buffer.get_chunksize();
    const auto refs = kkpair_size_t({ n_to_chunks, h_data.extent(0) });
    return Kokkos::subview(h_data, refs);
  }

 public:
  FSStoreArrayViaBuffer(FSStore& store, const std::vector<size_t> &chunkshape,
    const std::string_view name, const std::string_view units, const double scale_factor,
    const std::string_view dtype, const std::vector<std::string>& dims)
    : store(store), buffer(chunkshape), chunks(dims.size()), name(name), chunkshape(chunkshape) {
    /* number of names of dimensions must match number of chunks' dimensions,
    and chunksize according to buffer must match total size of a (shaped) chunk */
    assert(chunkshape.size() == dims.size());
    auto chunksize = size_t{1};
    for (const auto& c : chunkshape) { chunksize *= c; }
    assert(buffer.get_chunksize() == chunksize);

    /* make string of zarray metadata for array in zarr store (incomplete because missing shape) */
    const auto order = 'C';      // layout of bytes in each chunk of array in storage ('C' or 'F')
    const auto compressor = std::string{ "null" };   // compression of data when writing to store
    const auto fill_value = std::string{ "null" };   // fill value for empty datapoints in array
    const auto filters = std::string{ "null" };      // codec configurations for compression
    const auto zarr_format = '2';                    // storage spec. version 2

    partial_metadata = std::string(
      "\"chunks\": " +
      vec_to_string(chunkshape) +
      ",\n"
      "\"dtype\": \"" +
      std::string(dtype) +  // dtype = datatype stored in arrays e.g. "<f8"
      "\",\n"
      "\"order\": \"" +
      order +
      "\",\n"
      "\"compressor\": " +
      compressor +
      ",\n"
      "\"fill_value\": " +
      fill_value +
      ",\n"
      "\"filters\": " +
      filters +
      ",\n"
      "\"zarr_format\": " +
      zarr_format);

    /* make string of zattrs attribute information for array in zarr store */
    const auto arrayattrs = std::string(
      "{\n"
      "\"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dims) +                // names of each dimension of array
      ",\n"
      "\"units\": " +
      "\"" + std::string(units) + "\"" +    // units of coordinate being stored
      ",\n"
      "\"scale_factor\": " +
      std::to_string(scale_factor) +        // scale_factor of data
      ",\n}");

    write_zattrs_json(store, name, arrayattrs);
    write_zarray_json(store, name, zarr_metadata());
  };

  ~FSStoreArrayViaBuffer() {
    // write buffer to chunk if it isn't empty
    if (buffer.get_fill() > 0) {
      chunks.write_chunk(store, name, buffer);
      write_zarray_json(store, name, zarr_metadata());
    }
  };

  void write_data_to_zarr_array(const dualview_type::t_host h_data) {
    std::cout << "writing data to buffer / output\n";

    std::cout << "buffer size: " << buffer.get_chunksize() << "\n";
    std::cout << "buffer space: " << buffer.get_space() << "\n";
    std::cout << "initial data to add: " << h_data.extent(0) << "\n";

    auto h_data_rem = buffer.copy_to_buffer(h_data);

    std::cout << "after copy to buffer: " << h_data_rem.extent(0) << "\n";
    std::cout << "buffer space: " << buffer.get_space() << "\n";

    h_data_rem = write_chunks_in_store(h_data_rem);

    std::cout << "after writing to chunks: " << h_data_rem.extent(0) << "\n";
    std::cout << "buffer space: " << buffer.get_space() << "\n";

    h_data_rem = buffer.copy_to_buffer(h_data_rem);

    assert((h_data_rem.extent(0) == 0) && "there is leftover data remaining after writing array");

    std::cout << "final remaining data: " << h_data_rem.extent(0) << "\n";
    std::cout << "buffer space: " << buffer.get_space() << "\n";
  };
};

#endif    // ROUGHPAPER_ZARR_FSSTORE_ARRAY_HPP_
