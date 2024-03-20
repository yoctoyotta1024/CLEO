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
 * Last Modified: Wednesday 20th March 2024
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
    const auto it = datasetdims.find(dimnames.at(aa));
    reduced_arrayshape.push_back(it->second);
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
  // TODO(CB) move aliases to aliases.hpp
  ZarrArray<Store, T> zarr;
  std::vector<std::string> dimnames;

  /* set the shape of the array and its dimensions
  given the dimensions in the dataset and order of dims in dimnames */
  void set_arrayshape(const std::unordered_map<std::string, size_t>& datasetdims) {
    auto arrayshape_from_dims = std::vector<size_t>({});

    for (auto& dname : dimnames) {
      const auto it = datasetdims.find(dname);
      arrayshape_from_dims.push_back(it->second);
    }

    zarr.set_arrayshape(arrayshape_from_dims);
  }

  /**
   * @brief Writes chunks of data from a kokkos view in host memory to the Zarr array in a store.
   *
   * Calls write_chunks_to_store to write whole chunks of data into store. Then updates the shape
   * of each of the dimensions of the array to be consistent with the accumulated change in shape
   * of the array (due to the chunks that have been written). Note however, this function does not
   * (re-)write the .zarray json file's metadata for the shape of the array.
   * Function returns a (sub)view of the remaining data not written to a chunk (number of elements
   * in subview < chunksize).
   *
   * @param h_data Kokkos view of the data to write to the store in host memory.
   * @return The remaining data that was not written to chunks.
   */
  subviewh_buffer write_chunks_with_xarray_metadata(
      const std::unordered_map<std::string, size_t>& datasetdims, const subviewh_buffer h_data) {
    const auto shape_increment = write_chunks_to_store(h_data);

    if (shape_increment) {
      update_arraydims(datasetdims, store, shape_increment);
    }

    const auto n_to_chunks = nchunks_data * buffer.get_chunksize();
    const auto refs = kkpair_size_t({n_to_chunks, h_data.extent(0)});
    return Kokkos::subview(h_data, refs);
  }

 public:
  XarrayZarrArray(Store& store, const std::unordered_map<std::string, size_t>& datasetdims,
                  const std::string_view name, const std::string_view units,
                  const std::string_view dtype, const double scale_factor,
                  const std::vector<size_t>& chunkshape, const std::vector<std::string>& dimnames)
      : ZarrArray(store, name, dtype, chunkshape,
                  reduced_arrayshape_from_dims(datasetdims, dimnames)),
        dimnames(dimnames) {
    assert((chunkshape.size() == dimnames.size()) &&
           "number of named dimensions of array must match number dimensions of chunks");

    set_arrayshape(datasetdims);

    write_zattrs_json(store, name, make_xarray_metadata(units, scale_factor, dimnames));
  }

  /**
   * @brief Writes data from Kokkos view in host memory to chunks of a Zarr array in a store
   * via a buffer.
   *
   * Copies some data from the view to a buffer (until number of elements in buffer = chunksize),
   * then may write chunks of the array alongside the necessary metadata for a Zarr array into
   * a store. Finall copies any leftover data, number of elements < chunksize, into the buffer.
   * Assertion checks there is no remainng data unattended to.
   *
   * @param h_data The data in a Kokkos view in host memory which should be written to the array
   * in a store.
   */
  std::unordered_map<std::string, size_t> write_to_array(
      const std::unordered_map<std::string, size_t>& datasetdims, const viewh_buffer h_data) {
    auto h_data_rem = buffer.copy_to_buffer(h_data);

    h_data_rem = write_chunks_with_xarray_metadata(datasetdims, h_data_rem);

    h_data_rem = buffer.copy_to_buffer(h_data_rem);

    /* TODO(CB) docstrings */
    /* TODO(CB) update .zarray json file for Zarr metadata about shape of array according to
     * shape of dimensions as required for xarray and NetCDF compatibility -> eventually call
     something like "set_arrayshape(shape) and write_zarray_json(store, name, zarr_metadata());"
   */

    assert((h_data_rem.extent(0) == 0) && "there is leftover data remaining after writing array");

    return arraydims;
  };
};

#endif  // ROUGHPAPER_ZARR_XARRAY_ZARR_ARRAY_HPP_
