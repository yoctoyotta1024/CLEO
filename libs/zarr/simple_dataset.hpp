/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: simple_dataset.hpp
 * Project: zarr
 * Created Date: Monday 18th March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Structure to create a ZarrGroup which is xarray and netCDF compatible.
 */

#ifndef LIBS_ZARR_SIMPLE_DATASET_HPP_
#define LIBS_ZARR_SIMPLE_DATASET_HPP_

#include <Kokkos_Core.hpp>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "zarr/xarray_zarr_array.hpp"
#include "zarr/zarr_group.hpp"

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
class SimpleDataset {
 private:
  ZarrGroup<Store> group; /**< Reference to the zarr group object. */
  std::unordered_map<std::string, size_t>
      datasetdims; /**< map from name of each dimension in dataset to their size */

  /**
   * @brief Adds a dimension to the dataset.
   *
   * @param dim A pair containing the name and size of the dimension to be added.
   */
  void add_dimension(const std::pair<std::string, size_t> &dim) {
    datasetdims.insert({dim.first, dim.second});
  }

 public:
  /**
   * @brief Constructs a Dataset with the specified store object.
   *
   * This constructor initializes a Dataset with the provided store object by initialising a
   * ZarrGroup and writing some additional metatdata for Xarray and NetCDF.
   *
   * @param store The store object associated with the Dataset.
   */
  explicit SimpleDataset(Store &store) : group(store), datasetdims() {
    store[".zattrs"] =
        "{\n"
        "  \"creator\": \"Clara Bayley\",\n"
        "  \"title\": \"Dataset from CLEO is Xarray and NetCDF compatible Zarr Group of Arrays\""
        "\n}";
  }

  /**
   * @brief Returns the size of an existing dimension in the dataset.
   *
   * @param dimname A string for the name of the dimension in the dataset.
   * @return The size of (i.e. number of elements along) the dimension.
   */
  size_t get_dimension(const std::string &dimname) const { return datasetdims.at(dimname); }

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
   * @param scale_factor The scale factor of array data.
   * @param chunkshape The shape of the chunks of the array.
   * @param dimnames The names of each dimension of the array.
   * @return An instance of XarrayZarrArray representing the newly created array.
   */
  template <typename T>
  XarrayZarrArray<Store, T> create_array(const std::string_view name, const std::string_view units,
                                         const double scale_factor,
                                         const std::vector<size_t> &chunkshape,
                                         const std::vector<std::string> &dimnames) const {
    return XarrayZarrArray<Store, T>(group.store, datasetdims, name, units, scale_factor,
                                     chunkshape, dimnames);
  }

  /**
   * @brief Creates a new 1-D array for a coordinate of the dataset.
   *
   * @tparam T The data type of the coordinate array.
   * @param name The name of the new coordinate.
   * @param units The units of the coordinate.
   * @param scale_factor The scale factor of the coordinate data.
   * @param chunksize The size of each 1-D chunk of the coordinate array.
   * @param dimsize The initial size of the coordinate (number of elements along array).
   * @return An instance of XarrayZarrArray representing the newly created coordinate array.
   */
  template <typename T>
  XarrayZarrArray<Store, T> create_coordinate_array(const std::string_view name,
                                                    const std::string_view units,
                                                    const double scale_factor,
                                                    const size_t chunksize, const size_t dimsize) {
    add_dimension(std::pair<std::string, size_t>{name, dimsize});
    return create_array<T>(name, units, scale_factor, std::vector<size_t>{chunksize},
                           std::vector<std::string>{std::string(name)});
  }

  /**
   * @brief Creates a new ragged array in the dataset.
   *
   * @tparam T The data type of the array.
   * @param name The name of the new array.
   * @param units The units of the array data.
   * @param scale_factor The scale factor of array data.
   * @param chunkshape The shape of the chunks of the array.
   * @param dimnames The names of each dimension of the array.
   * @param sampledimname The names of the sample dimension of the array.
   * @return An instance of XarrayZarrArray representing the newly created ragged array.
   */
  template <typename T>
  XarrayZarrArray<Store, T> create_ragged_array(const std::string_view name,
                                                const std::string_view units,
                                                const double scale_factor,
                                                const std::vector<size_t> &chunkshape,
                                                const std::vector<std::string> &dimnames,
                                                const std::string_view sampledimname) const {
    return XarrayZarrArray<Store, T>(group.store, datasetdims, name, units, scale_factor,
                                     chunkshape, dimnames, sampledimname);
  }

  /**
   * @brief Creates a new raggedcount array in the dataset.
   *
   * @tparam T The data type of the array.
   * @param name The name of the new array.
   * @param units The units of the array data.
   * @param scale_factor The scale factor of array data.
   * @param chunkshape The shape of the chunks of the array.
   * @param dimnames The names of each dimension of the array.
   * @param sampledimname The names of the sample dimension of the array.
   * @return An instance of XarrayZarrArray representing the newly created raggedcount array.
   */
  template <typename T>
  XarrayZarrArray<Store, T> create_raggedcount_array(const std::string_view name,
                                                     const std::string_view units,
                                                     const double scale_factor,
                                                     const std::vector<size_t> &chunkshape,
                                                     const std::vector<std::string> &dimnames,
                                                     const std::string_view sampledimname) const {
    return XarrayZarrArray<Store, T>(group.store, datasetdims, name, units, scale_factor,
                                     chunkshape, dimnames, sampledimname);
  }

  /**
   * @brief Calls array's shape function to ensure the shape of the array matches
   * the dimensions of the dataset.
   *
   * @tparam T The data type of the array.
   * @param xzarr An instance of XarrayZarrArray representing the array.
   */
  template <typename T>
  void write_arrayshape(XarrayZarrArray<Store, T> &xzarr) const {
    xzarr.write_arrayshape(datasetdims);
  }

  /**
   * @brief Calls array's shape function to ensure the shape of the array matches
   * the dimensions of the dataset.
   *
   * @tparam T The data type of the array.
   * @param xzarr_ptr A shared pointer to the instance of XarrayZarrArray representing the array.
   */
  template <typename T>
  void write_arrayshape(const std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr) const {
    xzarr_ptr->write_arrayshape(datasetdims);
  }

  /**
   * @brief Calls array's shape function to write the shape of the array for a ragged array.
   *
   * @tparam T The data type of the array.
   * @param xzarr An instance of XarrayZarrArray representing the array.
   */
  template <typename T>
  void write_ragged_arrayshape(XarrayZarrArray<Store, T> &xzarr) const {
    xzarr.write_ragged_arrayshape();
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
                      const typename Buffer<T>::viewh_buffer h_data) const {
    xzarr.write_to_array(h_data);
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
   * @param xzarr_ptr A shared pointer to the instance of XarrayZarrArray representing the array.
   * @param h_data The data to be written to the array.
   */
  template <typename T>
  void write_to_array(const std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr,
                      const typename Buffer<T>::viewh_buffer h_data) const {
    xzarr_ptr->write_to_array(h_data);
    xzarr_ptr->write_arrayshape(datasetdims);
  }

  /**
   * @brief Writes 1 data element to a Zarr array in the dataset and calls
   * function to ensure the shape of the array matches the dimensions of the dataset.
   *
   * Function writes 1 data element to an array in the dataset and updates the metadata for the
   * shape of the array to ensure the size of each dimension of the array is consistent with the
   * dimensions of the dataset.
   *
   * @tparam T The data type of the array.
   * @param xzarr_ptr A shared pointer to the instance of XarrayZarrArray representing the array.
   * @param data The data element to be written to the array.
   */
  template <typename T>
  void write_to_array(const std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr,
                      const T data) const {
    xzarr_ptr->write_to_array(data);
    xzarr_ptr->write_arrayshape(datasetdims);
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
  void write_to_ragged_array(XarrayZarrArray<Store, T> &xzarr,
                             const typename Buffer<T>::viewh_buffer h_data) const {
    xzarr.write_to_array(h_data);
    xzarr.write_ragged_arrayshape();
  }
};

#endif  // LIBS_ZARR_SIMPLE_DATASET_HPP_
