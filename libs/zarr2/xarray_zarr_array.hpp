/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: xarray_zarr_array.hpp
 * Project: zarr2
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 28th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Structure to create a group obeying the Zarr storage specification version 2
 * (https://zarr.readthedocs.io/en/stable/spec/v2.html) in a given memory store.
 */

#ifndef LIBS_ZARR2_XARRAY_ZARR_ARRAY_HPP_
#define LIBS_ZARR2_XARRAY_ZARR_ARRAY_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <algorithm>
#include <cassert>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "./zarr_array.hpp"

/**
 * @brief Write attributes string to a store under a .zattrs key.
 *
 * Write some data under .zattrs key in store for an array called 'name'. The key and attrs data
 * could be anything, but for example .zattrs could be a json file in a file system store
 * (see FSStore) for the extra metadata which must exist in order to make Xarray and netCDF
 * happy when opening a Zarr dataset, e.g. by naming the dimensions of the
 * "{\"_ARRAY_DIMENSIONS\": [\"dimension_name\"]}";.
 *
 * @tparam Store The type of the store object where the metadata will be written.
 * @param store The store object where the metadata will be written.
 * @param name The name under which the .zarray key will be stored in the store.
 * @param metadata The metadata to write for the .zarray key.
 */
template <typename Store>
inline void write_zattrs_json(Store& store, std::string_view name, std::string_view attrs) {
  store[std::string(name) + "/.zattrs"] = attrs;
}

/**
 * @brief Calculate the reduced array shape of an array given the name of its dimensions and the
 * dataset's dimensions.
 *
 * Given the dimensions of a dataset and the names of the dimensions of an array, this function
 * calculates the reduced array shape by extracting the sizes of the dimensions from the dataset
 * which correspond to the provided dimension names for all except for the outermost dimension of
 * the array.
 *
 * @param datasetdims An unordered map containing the dimensions of the dataset.
 * @param dimnames A vector containing the names of the dimensions of the array (ordered from
 * outermost->innermost).
 * @return A vector of size_t representing the reduced array shape.
 */
inline std::vector<size_t> reduced_arrayshape_from_dims(
    const std::unordered_map<std::string, size_t>& datasetdims,
    const std::vector<std::string>& dimnames) {
  auto reduced_arrayshape = std::vector<size_t>({});

  for (size_t aa = 1; aa < dimnames.size(); ++aa) {
    const auto dsize = datasetdims.at(dimnames.at(aa));  // number of elements along a dimension
    reduced_arrayshape.push_back(dsize);
  }

  return reduced_arrayshape;
}

/**
 * @brief Converts vector of strings, e.g. for names of dimensions,into a single list written
 * as a string.
 *
 * @param dims The vector of strings to be converted.
 * @return A string representing the converted list.
 */
inline std::string vecstr_to_string(const std::vector<std::string>& dims) {
  auto dims_str = std::string{"["};
  for (const auto& d : dims) {
    dims_str += "\"" + d + "\",";
  }
  dims_str.pop_back();  // delete last ","
  dims_str += "]";
  return dims_str;
}

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make zarr array
 * compatible with Xarray and NetCDF.
 *
 * @param units The units of the array's coordinates.
 * @param scale_factor The scale factor of data.
 * @param dimnames The names of each dimension of the array.
 * @return A string representing the metadata.
 */
inline std::string make_xarray_metadata(const std::string_view units, const double scale_factor,
                                        const std::vector<std::string>& dimnames) {
  const auto zattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dimnames) +  // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(units) + "\"" +  // units of coordinate being stored
      ",\n"
      "  \"scale_factor\": " +
      std::to_string(scale_factor) +  // scale_factor of data
      "\n}");

  return zattrs;
}

/**
 * @brief Zarr array with additional metadata and functions to constrain the shape of array to the
 * shape of its dimensions in a dataset in order to ensure Zarr array is compatibile with NetCDF
 * and Xarray conventions.
 *
 * @tparam Store The type of the store object where the array will be stored.
 * @tparam T The data type of the array.
 */
template <typename Store, typename T>
class XarrayZarrArray {
 private:
  using viewh_buffer = Buffer<T>::viewh_buffer;
  ZarrArray<Store, T> zarr;           ///< zarr array in store
  std::vector<std::string> dimnames;  ///< ordered list of names of each dimenion of array
  std::vector<size_t> arrayshape;     ///< current size of the array along each of its dimensions
  size_t last_totnchunks;             ///< Number of chunks of array since arrayshape last written

  /**
   * @brief Sets shape of array along each dimension to be the same size as each of its dimensions
   * according to the dataset. Returns boolean for whether shape has changed (true) or not (false).
   *
   * The order of the dimensions in the array's shape is the order of dimensions in dimnames
   * (outermost -> innermost). Setting the shape to be conistent with the size of the dataset's
   * dimensions makes zarr array also consistent with Xarray and NetCDF conventions. Boolean
   * returns true if the shape of the array along any of its dimensions has changed.
   *
   * @param datasetdims Dictionary like object for the dimensions of the dataset.
   * @return bool = true if arrayshape along any of its dimensions has changed, false otherwise.
   */
  bool set_arrayshape(const std::unordered_map<std::string, size_t>& datasetdims) {
    auto ischange = std::vector<int>(arrayshape.size(), 0);

    for (size_t aa = 0; aa < dimnames.size(); ++aa) {
      const auto dsize = datasetdims.at(dimnames.at(aa));
      ischange.at(aa) = dsize - arrayshape.at(aa);
      arrayshape.at(aa) = dsize;
    }

    return std::any_of(ischange.begin(), ischange.end(), [](bool b) { return b; });
  }

 public:
  /**
   * @brief Constructs a new XarrayZarrArray object.
   *
   * @param store The store where the array will be stored.
   * @param datasetdims Dictionary like object for the dimensions of the dataset.
   * @param name The name of the array.
   * @param units The units of the array data.
   * @param dtype The data type of the array.
   * @param scale_factor The scale factor of array data.
   * @param chunkshape The shape of the array chunks.
   * @param dimnames The names of each dimension of the array (in order outermost->innermost).
   */
  XarrayZarrArray(Store& store, const std::unordered_map<std::string, size_t>& datasetdims,
                  const std::string_view name, const std::string_view units,
                  const std::string_view dtype, const double scale_factor,
                  const std::vector<size_t>& chunkshape, const std::vector<std::string>& dimnames)
      : zarr(store, name, dtype, chunkshape, true,
             reduced_arrayshape_from_dims(datasetdims, dimnames)),
        dimnames(dimnames),
        arrayshape(dimnames.size(), 0),
        last_totnchunks(0) {
    assert((chunkshape.size() == dimnames.size()) &&
           "number of named dimensions of array must match number dimensions of chunks");

    write_arrayshape(datasetdims);

    write_zattrs_json(store, name, make_xarray_metadata(units, scale_factor, dimnames));
  }

  ~XarrayZarrArray() { zarr.write_arrayshape(arrayshape); }

  /**
   * @brief Returns the name and size of the dimensions of the array (unordered).
   *
   * @return An unordered map containing the current dimensions of the array.
   */
  std::unordered_map<std::string, size_t> get_arraydims() const {
    auto arraydims = std::unordered_map<std::string, size_t>();
    for (size_t aa = 0; aa < dimnames.size(); ++aa) {
      arraydims.insert({dimnames.at(aa), arrayshape.at(aa)});
    }

    return arraydims;
  }

  /**
   * @brief Writes data from Kokkos view in host memory to chunks of a Zarr array in a store
   * via a buffer. Function does *not* write metadata to zarray .json file.
   *
   * Calls ZarrArray's write_to_array function to write data from Kokkos view in host memory to
   * chunks of a Zarr array in a store.
   *
   * @param h_data The data in a Kokkos view in host memory which should be written to the array
   * in a store.
   */
  void write_to_array(const viewh_buffer h_data) { zarr.write_to_array(h_data); };

  /**
   * @brief Writes 1 data element to a Zarr array in a store.
   * Function does *not* write metadata to zarray .json file.
   *
   * Calls ZarrArray's write_to_array function to write data to a Zarr array in a store (in chunks
   * via a buffer).
   *
   * @param data The data element which should be written to the array in a store.
   */
  void write_to_array(const T data) { zarr.write_to_array(data); };

  /**
   * @brief Sets shape of array along each dimension to be the same size as each of its dimensions
   * according to the dataset.
   *
   * The order of the dimensions in the array's shape is the order of dimensions in dimnames
   * (outermost -> innermost). Setting the shape to be conistent with the size of the dataset's
   * dimensions makes zarr array also consistent with Xarray and NetCDF conventions. If chunks have
   * been written since last writing of the arrayshape, and the shape of the array has changed, then
   * function also overwrites the .zarray json file with metadata containing the new shape of the
   * array.
   *
   * @param datasetdims Dictionary like object for the dimensions of the dataset.
   */
  void write_arrayshape(const std::unordered_map<std::string, size_t>& datasetdims) {
    auto ischange = set_arrayshape(datasetdims);

    if (last_totnchunks != zarr.get_totnchunks() && ischange) {
      zarr.write_arrayshape(arrayshape);
      last_totnchunks = zarr.get_totnchunks();
    }
  }
};

#endif  // LIBS_ZARR2_XARRAY_ZARR_ARRAY_HPP_
