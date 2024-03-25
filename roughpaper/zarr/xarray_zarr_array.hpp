/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: xarray_zarr_array.hpp
 * Project: zarr
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 25th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Structure to create a group obeying the Zarr storage specification version 2
 * (https://zarr.readthedocs.io/en/stable/spec/v2.html) in a given memory store.
 */

#ifndef ROUGHPAPER_ZARR_XARRAY_ZARR_ARRAY_HPP_
#define ROUGHPAPER_ZARR_XARRAY_ZARR_ARRAY_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
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
 * (see FSStore) for the extra metadata which must exist in order to make xarray and netCDF
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

/* converts vector of strings, e.g. for names of dimensions, into a single list
written as a string */
inline std::string vecstr_to_string(const std::vector<std::string>& dims) {
  auto dims_str = std::string{"["};
  for (const auto& d : dims) {
    dims_str += "\"" + d + "\",";
  }
  dims_str.pop_back();  // delete last ","
  dims_str += "]";
  return dims_str;
}

/* make string of array attributes metadata for .zattrs json for making zarr array
compatible with xarray and NetCDF */
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

/** Zarr array with additional metadata and constraining of shape of array to shape of dimensions
 * in order to ensure Zarr array is compatibile with NetCDF and Xarray
 * @param units The units of the array's coordinates.
 * @param dimnames The names of each dimension of the array.
 */
template <typename Store, typename T>
class XarrayZarrArray {
 private:
  using viewh_buffer = Buffer<T>::viewh_buffer;  // TODO(CB) move aliases to aliases.hpp
  ZarrArray<Store, T> zarr;                      ///< zarr array in store
  std::vector<std::string> dimnames;  ///< ordered list of names of each dimenion of array

  /* set the shape of the array and its dimensions
  given the dimensions in the dataset and order of dims in dimnames */
  void write_arrayshape(const std::unordered_map<std::string, size_t>& datasetdims) {
    auto arrayshape = std::vector<size_t>({});

    for (auto& dim : dimnames) {
      const auto dsize = datasetdims.at(dim);
      arrayshape.push_back(dsize);
    }

    zarr.write_arrayshape(arrayshape);
  }

  std::unordered_map<std::string, size_t> get_arraydims() const {
    auto arraydims = std::unordered_map<std::string, size_t>();
    auto arrayshape = zarr.get_arrayshape();
    for (size_t aa = 0; aa < dimnames.size(); ++aa) {
      arraydims.insert({dimnames.at(aa), arrayshape.at(aa)});
    }

    return arraydims;
  }

 public:
  XarrayZarrArray(Store& store, const std::unordered_map<std::string, size_t>& datasetdims,
                  const std::string_view name, const std::string_view units,
                  const std::string_view dtype, const double scale_factor,
                  const std::vector<size_t>& chunkshape, const std::vector<std::string>& dimnames)
      : zarr(store, name, dtype, chunkshape, true,
             reduced_arrayshape_from_dims(datasetdims, dimnames)),
        dimnames(dimnames) {
    assert((chunkshape.size() == dimnames.size()) &&
           "number of named dimensions of array must match number dimensions of chunks");

    write_arrayshape(datasetdims);  // overwrite zarr array shape with xarray dataset dimensions

    write_zattrs_json(store, name, make_xarray_metadata(units, scale_factor, dimnames));
  }

  /**
   * @brief Writes data from Kokkos view in host memory to chunks of a Zarr array in a store
   * via a buffer such that Zarr array in compatible with NetCDF and Xarray.
   *
   * Calls ZarrArray's write_to_array function to write data from Kokkos view in host memory to
   * chunks of a Zarr array in a store. Then overwrites the arrayshape and corresponding metadata
   * to ensure the shape of the array is consistent with the dimensions of the dataset, as required
   * by Xarray and NetCDF.
   *
   * @param h_data The data in a Kokkos view in host memory which should be written to the array
   * in a store.
   */
  void write_to_xarray_zarr_array(const std::unordered_map<std::string, size_t>& datasetdims,
                                  const viewh_buffer h_data) {
    zarr.write_to_array(h_data);
    write_arrayshape(datasetdims);  // overwrite zarr array shape with xarray dataset dimensions

    // TODO(CB) call write_arrayshape(datasetdims); after this function call in dataset once
    // dimensions of dataset have been consolidated (see WIP in dataset)

    // TODO(CB) docstrings
  };
};

#endif  // ROUGHPAPER_ZARR_XARRAY_ZARR_ARRAY_HPP_
