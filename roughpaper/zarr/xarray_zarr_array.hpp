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
 * Last Modified: Tuesday 19th March 2024
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

#include "./zarr_array.hpp"

/** Zarr array with additional metadata and constraining of shape of array to shape of dimensions
 * in order to ensure Zarr array is compatibile with NetCDF and Xarray
 * @param units The units of the array's coordinates.
 * @param dims The names of each dimension of the array.
 */
template <typename Store, typename T>
class XarrayZarrArray {
 private:
  ZarrArray<Store, T> zarr;
  std::vector<std::string> dims;

 public:
  XarrayZarrArray(Store& store, const std::string_view name, const std::string_view units,
                  const double scale_factor, const std::string_view dtype,
                  const std::vector<std::string>& dims, const std::vector<size_t>& chunkshape,
                  const std::vector<size_t>& reduced_arrayshape = std::vector<size_t>({}))
      : ZarrArray(store, name, dtype, chunkshape, reduced_arrayshape), dims(dims) {
    assert((chunkshape.size() == dims.size()) &&
           "number of named dimensions of array must match number dimensinos of chunks");

    /* make string of zattrs attribute information for array in zarr store */
    const auto arrayattrs = std::string(
        "{\n"
        "  \"_ARRAY_DIMENSIONS\": " +
        vecstr_to_string(dims) +  // names of each dimension of array
        ",\n"
        "  \"units\": " +
        "\"" + std::string(units) + "\"" +  // units of coordinate being stored
        ",\n"
        "  \"scale_factor\": " +
        std::to_string(scale_factor) +  // scale_factor of data
        "\n}");

    write_zattrs_json(store, name, arrayattrs);
  }

  void
  write_to_xarray_zarr_array(const viewh_buffer h_data) {
    auto h_data_rem = zarr.buffer.copy_to_buffer(h_data);

    h_data_rem = write_chunks_to_store(h_data_rem);
    chunk_metadata();

    h_data_rem = buffer.copy_to_buffer(h_data_rem);

    assert((h_data_rem.extent(0) == 0) && "there is leftover data remaining after writing array");
  };
};

#endif  // ROUGHPAPER_ZARR_XARRAY_ZARR_ARRAY_HPP_
