/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: zarr_array.hpp
 * Project: zarr2
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 27th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Class to write data to an array in a Zarr storage specification version 2
 * (https://zarr.readthedocs.io/en/stable/spec/v2.html) in a given memory store.
 */

#ifndef LIBS_ZARR2_ZARR_ARRAY_HPP_
#define LIBS_ZARR2_ZARR_ARRAY_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "./buffer.hpp"
#include "./chunks.hpp"

/**
 * @brief Given maximum chunk size 'maxchunksize' and length of inner dimension of one chunk of
 * array 'ndim1', function returns the largest possible chunk shape that has the length of its
 * inner dimension = ndim1.
 *
 * ndim1 must also be <= maxchunksize and to ensure good chunking, ndim1 should itself be completely
 * divisible by the final length of the inner dimension of the 2-D array.
 *
 * @param maxchunksize The maximum chunk size (maximum number of elements in chunk).
 * @param ndim1 The length of (number of elements along) the inner dimension of one chunk.
 * @return std::vector<size_t> The largest possible 2-D chunk shape.
 */
inline std::vector<size_t> good2Dchunkshape(const size_t maxchunksize, const size_t ndim1) {
  const auto shape0 = size_t{std::floor(maxchunksize / ndim1)};
  return {shape0, ndim1};
}

/**
 * @brief Write metadata string to a store under a .zarray key.
 *
 * write metadata under .zarray key in store for an array called 'name'. The key and metadata
 * could be anything, but for example .zarray could be a json file in a file system store
 * (see FSStore) for the metadata which must exist in order to decode chunks of an array according
 * to Zarr storage specification version 2 (https://zarr.readthedocs.io/en/stable/spec/v2.html),
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
 * @brief Converts a vector of integers into a single list written as a string.
 *
 * Given vector of a type convertible to a string with values [a, b, c, ..., z], function returns
 * the string "[a, b, c, ..., z]" with elements separated by commas and enclosed in square brackets.
 * Function is useful for converting vectors representing the shape of chunks and arrays etc. into
 * a string format for metadata json files.
 *
 * @param vals The vector values of a type convertible to a string.
 * @return A string representation of the vector.
 */
inline std::string vec_to_string(const std::vector<size_t>& vals) {
  auto vals_str = std::string{"["};
  for (const auto& v : vals) {
    vals_str += std::to_string(v) + ", ";
  }
  vals_str.erase(vals_str.size() - 2);  // delete last ", "
  vals_str += "]";
  return vals_str;
}

/**
 * @brief Generates part of the metadata for a Zarr array .zarray json file.
 *
 * This function constructs a string containing all the compulsory metadata of a Zarr array for its
 * .zarray json file, excluding the array's shape.
 *
 * @param chunkshape The shape of individual data chunks along each dimension.
 * @param dtype The data type stored in the arrays (e.g., "<f8").
 * @return A string view containing the partial metadata for the Zarr array.
 */
inline std::string make_part_zarrmetadata(const std::vector<size_t>& chunkshape,
                                          const std::string_view dtype) {
  const auto chunkshape_str = vec_to_string(chunkshape);  // shape of each chunk of array
  const auto order = 'C';  // layout of bytes in each chunk of array in storage ('C' or 'F')
  const auto compressor = std::string{"null"};  // compression of data when writing to store
  const auto fill_value = std::string{"null"};  // fill value for empty datapoints in array
  const auto filters = std::string{"null"};     // codec configurations for compression
  const auto zarr_format = '2';                 // storage spec. version 2

  const auto part_zarrmetadata = std::string("  \"chunks\": " + chunkshape_str +
                                             ",\n"
                                             "  \"dtype\": \"" +
                                             std::string(dtype) +
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
  return part_zarrmetadata;
}

/**
 * @brief A template class representing a Zarr array.
 *
 * This class provides functionality to write an array to a specified store via a buffer according
 * to the Zarr storage specification version 2 (https://zarr.readthedocs.io/en/stable/spec/v2.html).
 *
 * @tparam Store The type of store where the array will be stored.
 * @tparam T The data type stored in the arrays.
 */
template <typename Store, typename T>
class ZarrArray {
 private:
  using viewh_buffer = Buffer<T>::viewh_buffer;
  using subviewh_buffer = Buffer<T>::subviewh_buffer;
  Store& store;                   ///< store in which to write Zarr array
  std::string_view name;          ///< Name of array to write in store.
  size_t totnchunks;              ///< Total number of chunks of array written to store.
  size_t totndata;                ///< Total number of elements of data in array written to store.
  Chunks chunks;                  ///< Method to write chunks of array in store.
  Buffer<T> buffer;               ///< Buffer to hold data before writing chunks to store.
  std::string part_zarrmetadata;  ///< Metadata required for zarr array excluding array's shape
  bool is_backend;                ///< true if zarr array is a backend of something else e.g. xarray

  /**
   * @brief Generates the compulsory metadata for the Zarr array .zarray json file.
   *
   * @return A string containing the metadata for the Zarr array.
   */
  std::string zarr_metadata(const std::vector<size_t>& arrayshape) const {
    const auto metadata = std::string(
        "{\n"
        "  \"shape\": " +
        vec_to_string(arrayshape) + ",\n" + std::string(part_zarrmetadata) + "\n}");

    return metadata;
  }

  /**
   * @brief Get the shape of the array based on the number of data elements and chunks written in
   * the store.
   *
   * This method assumes that writing of chunks always fills inner dimensions first.
   * The array shape returned is always at least large enough to display all the elements of data
   * in the array so far along each dimension (i.e., arraysize >= totndata along each dimension of
   * the array).
   *
   * @return A vector representing the shape of the array.
   */
  std::vector<size_t> get_arrayshape() const {
    const auto chunkshape = chunks.get_chunkshape();
    const auto reducedarray_nchunks = chunks.get_reducedarray_nchunks();

    auto arrayshape = std::vector<size_t>(chunkshape.size(), 0);
    for (size_t aa = 1; aa < arrayshape.size(); ++aa) {
      const auto nchunks = vec_product(reducedarray_nchunks, aa);    // nchunks along inner dims
      const auto maxnchunks = (totnchunks + nchunks - 1) / nchunks;  // round up int division
      arrayshape.at(aa) = std::min(maxnchunks, reducedarray_nchunks.at(aa - 1)) * chunkshape.at(aa);
    }

    const auto reduced_arrayndata = size_t{std::max(vec_product(arrayshape, 1), size_t{1})};
    const auto wholeblocksize = (reduced_arrayndata * chunkshape.at(0));
    const auto whole_shape0 = (totndata / wholeblocksize) * chunkshape.at(0);

    const auto remainder_ndata = totndata - (whole_shape0 * reduced_arrayndata);
    const auto remainder_shape0 = std::min(remainder_ndata, chunkshape.at(0));
    arrayshape.at(0) = whole_shape0 + remainder_shape0;

    assert((totndata <= vec_product(arrayshape)) &&
           "elements of data must not be hiddden by array shape");
    return arrayshape;
  }

  /**
   * @brief Writes chunks of data from a kokkos view in host memory to the Zarr array in a store.
   *
   * First writes the buffer to a chunk of the array if it's full. Then writes whole chunks
   * directly from the Kokkos view if the view contains enough elements for whole chunk(s) to be
   * written. Then updates the shape of the array along its outermost dimension with the accumulated
   * change in shape of the array due to the chunks that have been written. Finally returns a
   * (sub)view of the remaining data not written to a chunk (number of elements in
   * subview < chunksize). Note that this function does not ensure the .zarray json file metadata
   * is kept up-to-date with changes to the arrayshape that may occur due to the increase in
   * number of elements of data written to the array during this function call.
   *
   * @param h_data Kokkos view of the data to write to the store in host memory.
   * @return The remaining data that was not written to chunks.
   */
  subviewh_buffer write_chunks_to_store(const subviewh_buffer h_data) {
    const auto csz = buffer.get_chunksize();

    if (buffer.get_space() == 0) {
      totnchunks = chunks.write_chunk<Store, T>(store, name, totnchunks, buffer);
    }

    const auto nchunks_data = size_t{h_data.extent(0) / csz};
    totnchunks = chunks.write_chunks<Store, T>(store, name, h_data, totnchunks, csz, nchunks_data);
    totndata = totnchunks * csz;

    const auto n_to_chunks = nchunks_data * csz;
    const auto refs = kkpair_size_t({n_to_chunks, h_data.extent(0)});
    return Kokkos::subview(h_data, refs);
  }

 public:
  /**
   * @brief Constructs a ZarrArray object.
   *
   * Initializes an empty Zarr array in the provided store in order to writes chunks of an array to
   * the store via a buffer. The assertions in this constructor ensure chunks have the same
   * number of dimensions for the array. The buffer is the size of exactly 1 chunk, and chunks'
   * shape is restricted such that the final array dimensions are exactly integer multiples of its
   * chunks along all but the outermost (0th) dimension of the array. Order of data written to
   * chunks is assumed to increment along innermost dimensions first.
   *
   * @param store The store where the array will be stored.
   * @param chunkshape The shape of individual data chunks along each dimension.
   * @param reduced_arrayshape The shape of the array along all but the outermost (0th) dimension.
   * @param name The name of the array.
   * @param scale_factor The scale factor of the data.
   * @param dtype The data type stored in the arrays (e.g., "<f8").
   */
  ZarrArray(Store& store, const std::string_view name, const std::string_view dtype,
            const std::vector<size_t>& chunkshape, const bool is_backend,
            const std::vector<size_t>& reduced_arrayshape = std::vector<size_t>({}))
      : store(store),
        name(name),
        totnchunks(0),
        totndata(0),
        chunks(chunkshape, reduced_arrayshape),
        buffer(vec_product(chunks.get_chunkshape())),
        part_zarrmetadata(make_part_zarrmetadata(chunkshape, dtype)),
        is_backend(is_backend) {
    assert((chunkshape.size() == reduced_arrayshape.size() + 1) &&
           "number of dimensions of chunks must match number of dimensions of array");

    /* Initial array shape is [0,0,0,...,0] (initially empty array along all dimensions) */
    write_arrayshape(get_arrayshape());
  };

  /**
   * @brief Destroys the ZarrArray object.
   *
   * Writes the buffer to a chunk of the array in the store if it isn't empty and issues a warning
   * if the data in buffer mismatches the array's expected dimensions. If the array is not a
   * backend (e.g. of an array in an xarray or NetCDF dataset), then the metadata for the
   * array's shape is also updated and warnings are issued if the array is incomplete.
   */
  ~ZarrArray() {
    if (buffer.get_fill() > 0) {
      if (buffer.get_fill() % vec_product(chunks.get_chunkshape(), 1) != 0) {
        std::cout << "WARNING: The number of data elements in the buffer is not completely "
                     "divisible by the number of elements in a chunk along its inner dimensions\n";
      }

      totndata = totnchunks * buffer.get_chunksize() + buffer.get_fill();
      totnchunks = chunks.write_chunk<Store, T>(store, name, totnchunks, buffer);
    }

    if (!(is_backend)) {
      write_arrayshape(get_arrayshape());

      const auto arrayshape = get_arrayshape();
      const auto reduced_arrayshape = chunks.get_reduced_arrayshape();
      for (size_t aa = 1; aa < arrayshape.size(); ++aa) {
        if (arrayshape.at(aa) < reduced_arrayshape.at(aa - 1)) {
          std::cout << "WARNING: array is not complete along inner dimension: " << aa << "\n";
        }
      }
      if (totndata < vec_product(arrayshape)) {
        std::cout << "WARNING: array is larger than total number of elements of data in it. Array"
                     "will have missing (i.e. null / nan) values.\n";
      }
    }
  }

  void write_arrayshape(const std::vector<size_t>& arrayshape) {
    assert((arrayshape.size() == chunks.get_chunkshape().size()) &&
           "number of dimensions of array must not change");
    write_zarray_json(store, name, zarr_metadata(arrayshape));
  }

  /**
   * @brief Writes data from Kokkos view in host memory to chunks of a Zarr array in a store
   * via a buffer and keep metadata in zarray .json file up-to-date with written chunks.
   *
   * First copies some data from the view to a buffer (until number of elements in
   * buffer = chunksize). Second writes any whole chunks of the array into a store. Thirdly updates
   * the .zarray json file for the Zarr metadata about the shape of the array accordingly. Finally
   * copies any leftover data, number of elements < chunksize, into the buffer.
   * Assertion checks there is no remainng data unattended to.
   *
   * @param h_data The data in a Kokkos view in host memory which should be written to the array in
   * a store.
   */
  void write_to_zarr_array(const viewh_buffer h_data) {
    auto h_data_rem = buffer.copy_to_buffer(h_data);

    h_data_rem = write_chunks_to_store(h_data_rem);
    write_arrayshape(get_arrayshape());  // ensure shape of array is up-to-date

    h_data_rem = buffer.copy_to_buffer(h_data_rem);

    assert((h_data_rem.extent(0) == 0) && "there is leftover data remaining after writing array");
  }

  /**
   * @brief Writes data from Kokkos view in host memory to chunks of Zarr array in a store
   * via a buffer. Function does *not* write metadata to zarray .json file.
   *
   * First copies some data from the view to a buffer (until number of elements in
   * buffer = chunksize), then writes any whole chunks of the array into a store. Finally
   * copies any leftover data, number of elements < chunksize, into the buffer.
   * Assertion checks there is no remainng data unattended to. Function useful when using zarr array
   * as backend of a dataset and/or you do not want to write metadata for the array when writing
   * data elements.
   *
   * @param h_data The data in a Kokkos view in host memory which should be written to the array in
   * a store.
   */
  void write_to_array(const viewh_buffer h_data) {
    auto h_data_rem = buffer.copy_to_buffer(h_data);

    h_data_rem = write_chunks_to_store(h_data_rem);

    h_data_rem = buffer.copy_to_buffer(h_data_rem);

    assert((h_data_rem.extent(0) == 0) && "there is leftover data remaining after writing array");
  }
};

#endif  // LIBS_ZARR2_ZARR_ARRAY_HPP_
