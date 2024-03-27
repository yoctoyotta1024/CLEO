/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: dataset.hpp
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
 * Structure to create a ZarrGroup which is xarray and netCDF compatible.
 */

#ifndef LIBS_ZARR2_DATASET_HPP_
#define LIBS_ZARR2_DATASET_HPP_

#include <Kokkos_Core.hpp>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "./xarray_zarr_array.hpp"
#include "./zarr_group.hpp"

/**
 * @brief A class representing a dataset made from a Zarr group (i.e. collection of Zarr arrays)
 * in a storage system.
 *
 * This class provides functionality to create a dataset as a group of arrays obeying the Zarr
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

  /**
   * @brief Adds a dimension to the dataset.
   *
   * @param dim A pair containing the name and size of the dimension to be added.
   */
  void add_dimension(const std::pair<std::string, size_t> &dim) {
    datasetdims.insert({dim.first, dim.second});
  }

  /**
   * @brief Sets the size of an existing dimension in the dataset.
   *
   * @param dim A pair containing the name of the dimension and its new size to be set.
   */
  void set_dimension(const std::pair<std::string, size_t> &dim) {
    datasetdims.at(dim.first) = dim.second;
  }

  /**
   * @brief Creates a new array in the dataset.
   *
   * @tparam T The data type of the array.
   * @param name The name of the new array.
   * @param units The units of the array data.
   * @param dtype The data type of the array.
   * @param scale_factor The scale factor of array data.
   * @param chunkshape The shape of the chunks of the array.
   * @param dimnames The names of each dimension of the array.
   * @return An instance of XarrayZarrArray representing the newly created array.
   */
  template <typename T>
  XarrayZarrArray<Store, T> create_array(const std::string_view name, const std::string_view units,
                                         const std::string_view dtype, const double scale_factor,
                                         const std::vector<size_t> &chunkshape,
                                         const std::vector<std::string> &dimnames) const {
    return XarrayZarrArray<Store, T>(group.store, datasetdims, name, units, dtype, scale_factor,
                                     chunkshape, dimnames);
  }

  /**
   * @brief Calls function to ensure the shape of the array matches the dimensions of the dataset.
   *
   * @tparam T The data type of the array.
   * @param xzarr An instance of XarrayZarrArray representing the array.
   */
  template <typename T>
  void write_arrayshape(XarrayZarrArray<Store, T> &xzarr) const {
    xzarr.write_arrayshape(datasetdims);
  }

  /**
   * @brief Writes data from Kokkos view in host memory to a Zarr array in the dataset and calls
   * function to ensure the shape of the array matches the dimensions of the dataset.
   *
   * Function writes data to an array in the dataset and updates the metadata for the shape of
   * the array to ensure the size of each dimension of the array is consistent with the
   * dimensions of the dataset.
   *
   * @tparam T The data type of the array.
   * @param xzarr An instance of XarrayZarrArray representing the array.
   * @param h_data The data to be written to the array.
   */
  template <typename T>
  void write_to_array(XarrayZarrArray<Store, T> &xzarr,
                      const Buffer<T>::viewh_buffer h_data) const {
    xzarr.write_to_array(h_data);
    xzarr.write_arrayshape(datasetdims);
  }
};

#endif  // LIBS_ZARR2_DATASET_HPP_
