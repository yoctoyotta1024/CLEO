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
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Structure to create a group obeying the Zarr storage specification version 2
 * (https://zarr.readthedocs.io/en/stable/spec/v2.html) in a given memory store.
 */

#ifndef LIBS_ZARR_XARRAY_ZARR_ARRAY_HPP_
#define LIBS_ZARR_XARRAY_ZARR_ARRAY_HPP_

#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <algorithm>
#include <cassert>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "configuration/communicator.hpp"
#include "zarr/xarray_metadata.hpp"
#include "zarr/zarr_array.hpp"

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
 * @param attrs The metadata to write for the .zarray key.
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
  ZarrArray<Store, T> zarr;          /**< zarr array in store */
  std::vector<std::string> dimnames; /**< ordered list of names of each dimenion of array */
  std::vector<size_t> arrayshape;    /**< current size of the array along each of its dimensions */
  size_t last_totnchunks;            /**< Number of chunks of array since arrayshape last written */
  MPI_Comm comm; /**< (YAC compatible) communicator for MPI domain decomposition */

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

  /**
   * @brief Sets shape of array along each dimension to be the same as shape according to zarr.
   *
   * Useful when writing a ragged arrray in a dataset (meaning length of dimensions if not length of
   * array)
   *
   */
  bool set_ragged_arrayshape() {
    const auto raggedarrayshape = std::vector<size_t>{zarr.get_totalndata()};
    const auto ischange = (arrayshape != raggedarrayshape);

    arrayshape = raggedarrayshape;

    return ischange;
  }

 public:
  /**
   * @brief Constructs a new XarrayZarrArray object.
   *
   * @param store The store where the array will be stored.
   * @param datasetdims Dictionary like object for the dimensions of the dataset.
   * @param name The name of the array.
   * @param units The units of the array data.
   * @param scale_factor The scale factor of array data.
   * @param chunkshape The shape of the array chunks.
   * @param dimnames The names of each dimension of the array (in order outermost->innermost).
   */
  XarrayZarrArray(Store& store, const std::unordered_map<std::string, size_t>& datasetdims,
                  const std::string_view name, const std::string_view units,
                  const double scale_factor, const std::vector<size_t>& chunkshape,
                  const std::vector<std::string>& dimnames)
      : zarr(store, name, chunkshape, true, reduced_arrayshape_from_dims(datasetdims, dimnames)),
        dimnames(dimnames),
        arrayshape(dimnames.size(), 0),
        last_totnchunks(0) {
    assert((chunkshape.size() == dimnames.size()) &&
           "number of named dimensions of array must match number dimensions of chunks");
    int my_rank;
    my_rank = init_communicator::get_comm_rank();

    if (my_rank == 0) {
      write_arrayshape(datasetdims);
      write_zattrs_json(store, name, xarray_metadata<T>(units, scale_factor, dimnames));
    }
  }

  /**
   * @brief Constructs a new XarrayZarrArray object with additional variable called
   * "sample_dimension" in the metadata .zattrs json and initially no set arrayshape.
   *
   * @param store The store where the array will be stored.
   * @param datasetdims Dictionary like object for the dimensions of the dataset.
   * @param name The name of the array.
   * @param units The units of the array data.
   * @param scale_factor The scale factor of array data.
   * @param chunkshape The shape of the array chunks.
   * @param dimnames The names of each dimension of the array (in order outermost->innermost).
   * @param sampledimname The name of the dimension the ragged count samples.
   */
  XarrayZarrArray(Store& store, const std::unordered_map<std::string, size_t>& datasetdims,
                  const std::string_view name, const std::string_view units,
                  const double scale_factor, const std::vector<size_t>& chunkshape,
                  const std::vector<std::string>& dimnames, const std::string_view sampledimname)
      : zarr(store, name, chunkshape, true, reduced_arrayshape_from_dims(datasetdims, dimnames)),
        dimnames(dimnames),
        arrayshape(dimnames.size(), 0),
        last_totnchunks(0) {
    assert((chunkshape.size() == dimnames.size()) &&
           "number of named dimensions of array must match number dimensions of chunks");
    int my_rank;
    my_rank = init_communicator::get_comm_rank();
    if (my_rank == 0) {
      write_zattrs_json(store, name,
                        xarray_metadata<T>(units, scale_factor, dimnames, sampledimname));
    }
  }

  ~XarrayZarrArray() {
    int my_rank;
    my_rank = init_communicator::get_comm_rank();
    if (my_rank == 0) zarr.write_arrayshape(arrayshape);
  }

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

  std::vector<std::string> get_dimnames() const { return dimnames; }

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

  /**
   * @brief Sets shape of array along each dimension to be the same a expected for a 1-D ragged
   * array.
   *
   * Expected shape is 1-D array with size of the total number of elements written to a zarr array.
   *
   * If chunks have been written since last writing of the arrayshape, then function also overwrites
   * the .zarray json file with metadata containing the new shape of the array.
   *
   */
  void write_ragged_arrayshape() {
    auto ischange = set_ragged_arrayshape();

    if (last_totnchunks != zarr.get_totnchunks() && ischange) {
      zarr.write_arrayshape(arrayshape);
      last_totnchunks = zarr.get_totnchunks();
    }
  }
};

#endif  // LIBS_ZARR_XARRAY_ZARR_ARRAY_HPP_
