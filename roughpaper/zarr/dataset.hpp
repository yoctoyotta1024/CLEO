/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: dataset.hpp
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
 * Structure to create a ZarrGroup which is xarray and netCDF compatible.
 */

#ifndef ROUGHPAPER_ZARR_DATASET_HPP_
#define ROUGHPAPER_ZARR_DATASET_HPP_

#include <Kokkos_Core.hpp>

#include "./xarray_zarr_array.hpp"
#include "./zarr_group.hpp"

/**
 * @brief A class representing a dataset made from a Zarr group (i.e. collection of Zarr arrays)
 * in a storage system.
 *
 * This class provides functionality to create a dataset  as a group of arrays obeying the Zarr
 * storage specification version 2 (https://zarr.readthedocs.io/en/stable/spec/v2.html) that is also
 * compatible with Xarray and NetCDF.
 *
 * @tparam Store The type of the store object used by the dataset.
 */
template <typename Store>
class Dataset {
 private:
  ZarrGroup<Store> group;  ///< Reference to the zarr group object.
  std::unordered_map<std::string, size_t>
      datasetdims;  ///< map from name of each dimension in dataset to their size

 public:
  /**
   * @brief Constructs a Dataset with the specified store object.
   *
   * This constructor initializes a Dataset with the provided store object by initialising a
   * ZarrGroup and writing some additional metatdata for Xarray and NetCDF.
   *
   * @param store The store object associated with the Dataset.
   */
  explicit Dataset(Store &store) : group(store), datasetdims() {
    store[".zattrs"] =
        "{\n"
        "  \"creator\": \"Clara Bayley\",\n"
        "  \"title\": \"Dataset from CLEO is Xarray and NetCDF compatible Zarr Group of Arrays\""
        "\n}";
  }

  void add_dimension(const std::pair<std::string, std::size_t> &dim) {
    datasetdims.insert({dim.first, dim.second});
  }

  void set_dimension(const std::pair<std::string, std::size_t> &dim) {
    datasetdims.at(dim.first) = dim.second;
  }

  template <typename T>
  XarrayZarrArray<Store, T> create_array(const std::string_view name, const std::string_view units,
                                         const std::string_view dtype, const double scale_factor,
                                         const std::vector<size_t> &chunkshape,
                                         const std::vector<std::string> &dimnames) const {
    return XarrayZarrArray<Store, T>(group.store, datasetdims, name, units, dtype, scale_factor,
                                     chunkshape, dimnames);
  }

  template <typename T>
  void write_to_array(XarrayZarrArray<Store, T> &xzarr,
                      const Buffer<T>::viewh_buffer h_data) const {
    xzarr.write_to_xarray_zarr_array(datasetdims, h_data);
  }
};

#endif  // ROUGHPAPER_ZARR_DATASET_HPP_
