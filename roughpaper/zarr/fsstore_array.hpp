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
 * Last Modified: Wednesday 13th March 2024
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
 public:
  const size_t chunksize;

  explicit Buffer(const size_t chunksize) : chunksize(chunksize), fill(0),
    buffer("buffer", chunksize) {
    reset_buffer();
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

 private:
  size_t fill;
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
};


struct ArrayChunks {
 private:
  const std::vector<size_t> chunkshape;    // shape of chunks of array along each dimension

  /* converts vector of integers for chunkcount into string to use to name a chunk */
  std::string chunkcount_to_string() {
    auto chunk_str = std::string{ "" };
    for (const auto& c : chunkcount) { chunk_str += std::to_string(c) + "."; }
    chunk_str.pop_back();   // delete last "."

    return chunk_str;
  }

  void update_chunks() {
    for (size_t aa = 0; aa < chunkshape.size(); ++aa) {
      chunkcount.at(0) += 1;             //+ ".0"; TODO(CB) deal with multi-D chunks
      shape.at(aa) += chunkshape.at(aa);
    }
  }

 public:
  std::vector<size_t> chunkcount;          // number of chunks along each dimension written in store
  std::vector<size_t> shape;               // number of elements in array along each dimension

  explicit ArrayChunks(const std::vector<size_t>& cs) : chunkshape(cs),
    chunkcount(std::vector<size_t>(cs.size(), 0)),
    shape(std::vector<size_t>(cs.size(), 0)) {}

  /* returns chunksize, the total number of elements in a chunk, which is the product of the number
  of elements along each dimension ie. the product of a chunk's shape. */
  size_t get_chunksize() {
    auto chunksize = size_t{1};
    for (const auto& c : chunkshape) { chunksize *= c; }
    return chunksize;
  }

  std::vector<size_t> get_chunkshape() {
    return chunkshape;
  }

  void write_chunk(FSStore& store, std::string_view name, Buffer buffer) {
    const auto chunk_str = chunkcount_to_string();
    buffer.write_buffer_to_chunk(store, name, chunk_str);
    update_chunks();
  }

  void write_chunk(FSStore& store, std::string_view name, const subview_type h_data_chunk) {
    const auto chunk_str = chunkcount_to_string();
    std::cout << "--> writing h_data to chunk: " << chunk_str << "\n";
    // write_data_to_chunk(store, name, chunk_str, h_data_chunk);   // TODO(CB) write subview chunk
    update_chunks();
  }
};

class FSStoreArrayViaBuffer {
 private:
  FSStore& store;                 // file system store satisfying zarr store specificaiton v2
  ArrayChunks chunks;             // information about chunks written in FSStore array
  Buffer buffer;                  // buffer for holding data before writing chunks to FSStore array
  std::string_view name;          // name to call variable being stored
  std::string partial_metadata;   // metadata excluding shape required for zarr array

  /* make string of metadata for array in zarr store */
  std::string zarr_metadata() {
    const auto metadata = std::string(
      "{\n"
      "\"shape\": " +
      vec_to_string(chunks.shape) +
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
    const auto nchunks_data = size_t{ h_data.extent(0) / buffer.chunksize };
    std::cout << "nchunks from h_data: " << nchunks_data << "\n";
    for (size_t bb = 0; bb < nchunks_data; ++bb) {
      const auto start = size_t{ bb * buffer.chunksize };
      const auto end = size_t{start + buffer.chunksize};
      const auto refs = kkpair_size_t({ start, end });
      chunks.write_chunk(store, name, Kokkos::subview(h_data, refs));
    }

    // update zarry json with new metadata
    write_zarray_json(store, name, zarr_metadata());

    // return remainder of data not written to chunks
    const auto n_to_chunks = nchunks_data * buffer.chunksize;
    const auto refs = kkpair_size_t({ n_to_chunks, h_data.extent(0) });
    return Kokkos::subview(h_data, refs);
  }

 public:
  FSStoreArrayViaBuffer(FSStore& store, const std::vector<size_t> &chunkshape,
    const std::string_view name, const std::string_view units, const double scale_factor,
    const std::string_view dtype, const std::vector<std::string>& dims)
    : store(store), chunks(chunkshape), buffer(chunks.get_chunksize()), name(name) {
    /* number of names of dimensions must match number of chunks' dimensions,
    and size of buffer must match size of an array chunk */
    assert(dims.size() == chunkshape.size());
    assert(buffer.chunksize == chunks.get_chunksize());

    /* make string of zarray metadata for array in zarr store (incomplete because missing shape) */
    const auto order = 'C';      // layout of bytes in each chunk of array in storage ('C' or 'F')
    const auto compressor = std::string{ "null" };   // compression of data when writing to store
    const auto fill_value = std::string{ "null" };   // fill value for empty datapoints in array
    const auto filters = std::string{ "null" };      // codec configurations for compression
    const auto zarr_format = '2';                    // storage spec. version 2

    partial_metadata = std::string(
      "\"chunks\": " +
      vec_to_string(chunks.get_chunkshape()) +
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
    if (buffer.get_space() < buffer.chunksize) {
      chunks.write_chunk(store, name, buffer);
    }
  };

  void write_data_to_zarr_array(const dualview_type::t_host h_data) {
    std::cout << "writing data to buffer / output\n";

    std::cout << "buffer size: " << buffer.chunksize << "\n";
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
