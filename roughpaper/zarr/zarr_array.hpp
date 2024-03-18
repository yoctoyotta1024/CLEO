/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: zarr_array.hpp
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
 * Class to write data to an array in a Zarr storage specification version 2
 * (https://zarr.readthedocs.io/en/stable/spec/v2.html) in a given memory store.
 */


#ifndef ROUGHPAPER_ZARR_ZARR_ARRAY_HPP_
#define ROUGHPAPER_ZARR_ZARR_ARRAY_HPP_

#include "./buffer.hpp"
#include "./chunks.hpp"

template <typename Store, typename T>
class ZarrArray {
 private:
  // TODO(CB) (1st move then) use aliases in aliases.hpp
  using viewh_buffer = Buffer<T>::viewh_buffer;
  using subviewh_buffer = Buffer<T>::subviewh_buffer;
  std::string_view name             ///< Name of array to write in store.
  size_t totnchunks                 ///< Total number of chunks of array written to store.
  std::vector<size_t> arrayshape;   ///< Number of elements in array along each dimension in store.
  Chunks chunks;                    ///< Method to write chunks of array in store.
  Buffer<T> buffer;                 ///< Buffer to hold data before writing chunks to store.
  std::string partial_metadata;     ///< Metadata required for zarr array excluding array's shape

  /* make string of metadata for array in zarr store */
  std::string zarr_metadata() const {
    const auto metadata = std::string(
      "{\n"
      "  \"shape\": " +
      vec_to_string(arrayshape) +
      ",\n" +
      std::string(partial_metadata) +
      "\n}");

    return metadata;
  }

  /* boolean returns increment to add to first dimension of array shape if chunk number indicates
  shape of array (and metadata) should be updated because the reduced array shape is complete,
  i.e. because the array along all but its outermost dimension is full of data elements */
  size_t arrayshape_change(const size_t chunk_num, const size_t shape_increment) const {
    if (chunk_num % vec_product(chunks.get_reducedarray_nchunks()) == 0) {
      return shape_increment;   // true
    } else {
      return 0;   // false
    }
  }

  /* increment shape of outermost dimension of N-Dimensional array and update the array's metadata
  .zarray json correspondingly. Function should only be called if is_arrayshape_change flag is
  curerntly true */
  void update_arrayshape(Store& store, const size_t shape_increment) {
    arrayshape.at(0) += shape_increment;   // increase in shape of outermost dimension
    write_zarray_json(store, name, zarr_metadata());   // update metadata
  }

  subviewh_buffer write_chunks_to_store(Store& store, const subviewh_buffer h_data) {
    auto shape_increment = size_t{ 0 };

    // write buffer to chunk if it's full
    if (buffer.get_space() == 0) {
      shape_increment += arrayshape_change(totnchunks, chunks.get_chunkshape().at(0));
      totnchunks = chunks.write_chunk<T>(store, name, totnchunks, buffer);
    }

    // write whole chunks of h_data_remaining
    const auto nchunks_data = size_t{ h_data.extent(0) / buffer.get_chunksize() };
    for (size_t nn = 0; nn < nchunks_data; ++nn) {
      const auto csz = buffer.get_chunksize();
      const auto refs = kkpair_size_t({ nn * csz, (nn + 1) * csz });
      shape_increment += arrayshape_change(totnchunks, chunks.get_chunkshape().at(0));
      totnchunks = chunks.write_chunk<T>(store, name, totnchunks, Kokkos::subview(h_data, refs));
    }

    // update shape of array if neccessary
    if (shape_increment) {
      update_arrayshape(store, shape_increment);
    }

    // return remainder of data not written to chunks
    const auto n_to_chunks = nchunks_data * buffer.get_chunksize();
    const auto refs = kkpair_size_t({ n_to_chunks, h_data.extent(0) });
    return Kokkos::subview(h_data, refs);
  }

 public:
  /**
  * @brief Writes a zarr array to a specified store via a buffer.
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
  ZarrArray(Store& store, const std::string_view name, const std::string_view units,
    const double scale_factor, const std::string_view dtype, const std::vector<std::string>& dims,
    const std::vector<size_t>& chunkshape,
    const std::vector<size_t>& reduced_arrayshape = std::vector<size_t>({}))
    : name(name), totnchunks(0), arrayshape(chunkshape.size(), 0),
    chunks(chunkshape, reduced_arrayshape), buffer(vec_product(chunks.get_chunkshape())) {
    /* number of names of dimensions must match number of dimensions of chunks */
    assert((chunkshape.size() == dims.size()) &&
      "number of named dimensions of array must match number dimensinos of chunks");

    /* number of dimensions for number of chunks must match number of dimensions of array */
    assert((chunkshape.size() == arrayshape.size()) &&
      "number of dimensions of chunks must match number of dimensions of array");

    /* make string of zarray metadata for array in zarr store (incomplete because missing shape) */
    const auto order = 'C';        // layout of bytes in each chunk of array in storage ('C' or 'F')
    const auto compressor = std::string{ "null" };   // compression of data when writing to store
    const auto fill_value = std::string{ "null" };   // fill value for empty datapoints in array
    const auto filters = std::string{ "null" };      // codec configurations for compression
    const auto zarr_format = '2';                    // storage spec. version 2


    /* set array shape along all but the outermost dimension to the number of
    elements given by the reduced array shape along that same dimension */
    for (size_t aa = 0; aa < reduced_arrayshape.size(); ++aa) {
      arrayshape.at(aa + 1) = reduced_arrayshape.at(aa);
    }

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

    write_zarray_json(store, name, zarr_metadata());

    // TODO(CB) move xarray metadata to dataset
    // /* make string of zattrs attribute information for array in zarr store */
    // const auto arrayattrs = std::string(
    //   "{\n"
    //   "  \"_ARRAY_DIMENSIONS\": " +
    //   vecstr_to_string(dims) +                // names of each dimension of array
    //   ",\n"
    //   "  \"units\": " +
    //   "\"" + std::string(units) + "\"" +    // units of coordinate being stored
    //   ",\n"
    //   "  \"scale_factor\": " +
    //   std::to_string(scale_factor) +        // scale_factor of data
    //   "\n}");

    // write_zattrs_json(store, name, arrayattrs);
  };

  ~ZarrArray() {
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

      auto shape_increment = buffer.get_fill() / reduced_arraysize;
      shape_increment = arrayshape_change(totnchunks, shape_increment);
      totnchunks = chunks.write_chunk<T>(store, name, totnchunks, buffer);
      update_arrayshape(store, shape_increment);
    }
  };

  void write_data_to_zarr_array(const viewh_buffer h_data) {
    auto h_data_rem = buffer.copy_to_buffer(h_data);

    h_data_rem = write_chunks_to_store(h_data_rem);

    h_data_rem = buffer.copy_to_buffer(h_data_rem);

    assert((h_data_rem.extent(0) == 0) && "there is leftover data remaining after writing array");
  };
};

#endif   // ROUGHPAPER_ZARR_ZARR_ARRAY_HPP_
