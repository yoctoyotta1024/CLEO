/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: collective_dataset.hpp
 * Project: zarr
 * Created Date: Friday 02nd August 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors: Wilton Jaciel Loch
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Distributed memory enabled structure to create a ZarrGroup which is xarray
 * and netCDF compatible. All processes have data, only one process creates the
 * output arrays, receives data, organizes it, and writes to the filesystem.
 */

#ifndef LIBS_ZARR_COLLECTIVE_DATASET_HPP_
#define LIBS_ZARR_COLLECTIVE_DATASET_HPP_

#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "configuration/communicator.hpp"
#include "zarr/xarray_zarr_array.hpp"
#include "zarr/zarr_group.hpp"

/**
 * @brief A class representing a dataset made from a Zarr group (i.e. collection of Zarr arrays)
 * in a storage system.
 *
 * This class provides functionality to create a dataset as a group of arrays obeying the Zarr
 * storage specification version 2 (https://zarr.readthedocs.io/en/stable/spec/v2.html) that is
 * also compatible with Xarray and NetCDF.
 *
 * @tparam Store The type of the store object used by the dataset.
 */
template <typename Store, typename Decomposition>
class CollectiveDataset {
 private:
  /**< Reference to the zarr group object. */
  ZarrGroup<Store> group;
  /**< map from name of each dimension in dataset to their size */
  std::unordered_map<std::string, size_t> datasetdims;
  Decomposition decomposition;
  std::shared_ptr<std::vector<unsigned int>> global_superdroplet_ordering;

  /**< map from name of each dimension in dataset to their size */
  std::unordered_map<std::string, std::vector<size_t>> distributed_datasetdims;
  int my_rank, comm_size;
  MPI_Comm comm; /**< (YAC compatible) communicator for MPI domain decomposition */

  /**
   * @brief Collects the distributed process-local size of dimensions
   *
   * @param dim A pair with the dimension name and local size
   */
  void collect_distributed_dim_size(const std::pair<std::string, size_t> &dim) {
    my_rank = init_communicator::get_comm_rank();
    comm_size = init_communicator::get_comm_size();
    std::vector<size_t> distributed_sizes(comm_size);

    MPI_Gather(&(dim.second), 1, MPI_UNSIGNED_LONG, distributed_sizes.data(), 1, MPI_UNSIGNED_LONG,
               0, comm);

    // For process 0, save the distributed sizes of all processes to avoid extra
    // exchanges in each write
    if (my_rank == 0) {
      if (distributed_datasetdims.contains(dim.first))
        distributed_datasetdims.at(dim.first) = distributed_sizes;
      else
        distributed_datasetdims.insert({dim.first, distributed_sizes});
    }
  }

  /**
   * @brief Collects the distributed process-local data for a write
   *
   * @param data A Kokkos view containing the data to be collected
   * @param dimnames The names of the dimensions related to the array
   */
  template <typename T>
  Kokkos::View<T *, HostSpace> collect_global_data(Kokkos::View<T *, HostSpace> data,
                                                   std::vector<std::string> dimnames) const {
    size_t innermost_dimension = dimnames.size() - 1;
    if (my_rank == 0) {
      std::vector<int> receive_displacements(comm_size, 0), receive_counts(comm_size, 0);
      size_t global_size = datasetdims.at(dimnames[innermost_dimension]);
      Kokkos::View<T *, HostSpace> global_data("global_output_data", global_size);

      // Calculate the receive counts and displacements to receive the data
      for (int i = 0; i < comm_size; i++) {
        receive_counts[i] =
            static_cast<int>(distributed_datasetdims.at(dimnames[innermost_dimension])[i]);
        if (i > 0)
          receive_displacements[i] =
              receive_displacements[i - 1] +
              distributed_datasetdims.at(dimnames[innermost_dimension])[i - 1];
      }

      // If the exchanges are for gridbox data extra actions are needed
      if (dimnames[innermost_dimension] == "gbxindex") {
        // If there is only the gridbox dimension, then this is the gridbox
        // index array and can be trivially constructed without exchanges
        if (dimnames.size() == 1) {
          for (size_t i = 0; i < decomposition.get_total_global_gridboxes(); i++)
            global_data[i] = i;
        } else {
          // All other gridbox-defined data must be correctly ordered to follow
          // the global gridbox ordering
          std::vector<T> receive_target(global_size);
          collect_global_array(receive_target.data(), data.data(),
                               distributed_datasetdims.at(dimnames[innermost_dimension])[0],
                               receive_counts.data(), receive_displacements.data());
          correct_gridbox_data(dimnames[innermost_dimension], global_data.data(),
                               receive_target.data());
        }
      } else {
        collect_global_array(global_data.data(), data.data(),
                             distributed_datasetdims.at(dimnames[innermost_dimension])[0],
                             receive_counts.data(), receive_displacements.data());
      }

      return global_data;
    } else {
      // Don't perform an exchange for the gridbox index array
      if (dimnames.size() != 1)
        collect_global_array(nullptr, data.data(), datasetdims.at(dimnames[innermost_dimension]),
                             nullptr, nullptr);

      return Kokkos::View<T *, HostSpace>();
    }
  }

  /**
   * @brief Correcly orders global data following the global gridbox order
   *
   * Given a source and a target array, correctly orders data from the source on
   * the target following the global gridbox ordering. Should be called only on
   * process 0.
   *
   * @param dimension The dimension name which should be the gridboxes
   * @param source The array to take the data from
   * @param target The array to write the data to, according to the global
   * gridbox order
   */
  template <typename T>
  void correct_gridbox_data(std::string dimension, T *target, T *source) const {
    size_t offset = 0;
    for (int process = 0; process < comm_size; process++) {
      for (size_t i = 0; i < distributed_datasetdims.at(dimension)[process]; i++) {
        size_t global_gridbox_index = decomposition.local_to_global_gridbox_index(i, process);
        target[global_gridbox_index] = source[offset + i];
      }
      offset += distributed_datasetdims.at(dimension)[process];
    }
  }

  /**
   * @brief Wrapper for MPI gatherv call for a float array
   */
  void collect_global_array(float *target, float *local_source, int local_size, int *receive_counts,
                            int *receive_displacements) const {
    MPI_Gatherv(local_source, local_size, MPI_FLOAT, target, receive_counts, receive_displacements,
                MPI_FLOAT, 0, comm);
  }

  /**
   * @brief Wrapper for MPI gatherv call for a unsigned int array
   */
  void collect_global_array(unsigned int *target, unsigned int *local_source, int local_size,
                            int *receive_counts, int *receive_displacements) const {
    MPI_Gatherv(local_source, local_size, MPI_UNSIGNED, target, receive_counts,
                receive_displacements, MPI_UNSIGNED, 0, comm);
  }

  /**
   * @brief Wrapper for MPI gatherv call for a unsigned int array
   */
  void collect_global_array(size_t *target, size_t *local_source, int local_size,
                            int *receive_counts, int *receive_displacements) const {
    MPI_Gatherv(local_source, local_size, MPI_UNSIGNED_LONG, target, receive_counts,
                receive_displacements, MPI_UNSIGNED_LONG, 0, comm);
  }

  /**
   * @brief Adds a dimension to the dataset.
   *
   * @param dim A pair containing the name and size of the dimension to be added.
   */
  void add_dimension(const std::pair<std::string, size_t> &dim) {
    collect_distributed_dim_size(dim);
    size_t dim_size = dim.second;

    // The time dimension has the global size in all processes and therefore
    // should not be accumulated
    if (my_rank == 0 && dim.first != "time")
      dim_size = std::accumulate(distributed_datasetdims.at(dim.first).begin(),
                                 distributed_datasetdims.at(dim.first).end(), 0);

    datasetdims.insert({dim.first, dim_size});
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
  explicit CollectiveDataset(Store &store) : group(store), datasetdims() {
    store[".zattrs"] =
        "{\n"
        "  \"creator\": \"Clara Bayley\",\n"
        "  \"title\": \"Dataset from CLEO is Xarray and NetCDF compatible Zarr Group of Arrays\""
        "\n}";
    global_superdroplet_ordering = std::make_shared<std::vector<unsigned int>>();
    comm = init_communicator::get_communicator();
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
    collect_distributed_dim_size(dim);
    size_t dim_size = dim.second;

    // The time dimension has the global size in all processes and therefore
    // should not be accumulated
    if (my_rank == 0 && dim.first != "time")
      dim_size = std::accumulate(distributed_datasetdims.at(dim.first).begin(),
                                 distributed_datasetdims.at(dim.first).end(), 0);

    datasetdims.at(dim.first) = dim_size;
  }

  /**
   * @brief Sets the decomposition maps for correctly writing data out
   *
   * @param decomposition A Decomposition instance for CLEO's domain decomposition
   */
  void set_decomposition(Decomposition decomposition) { this->decomposition = decomposition; }

  /**
   * @brief Sets the maximum number of superdroplets for data allocation, comes from the config file
   *
   * @param max_superdroplets The maximum number of superdroplets of the model
   */
  void set_max_superdroplets(unsigned int max_superdroplets) {
    constexpr unsigned int fill_value = std::numeric_limits<unsigned int>::max();
    global_superdroplet_ordering.get()->resize(max_superdroplets, fill_value);
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
    if (my_rank == 0) xzarr.write_arrayshape(datasetdims);
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
    if (my_rank == 0) xzarr_ptr->write_arrayshape(datasetdims);
  }

  /**
   * @brief Calls array's shape function to write the shape of the array for a ragged array.
   *
   * @tparam T The data type of the array.
   * @param xzarr An instance of XarrayZarrArray representing the array.
   */
  template <typename T>
  void write_ragged_arrayshape(XarrayZarrArray<Store, T> &xzarr) const {
    if (my_rank == 0) xzarr.write_ragged_arrayshape();
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
    auto global_data = collect_global_data(h_data, xzarr.get_dimnames());
    if (my_rank == 0) {
      xzarr.write_to_array(global_data);
      xzarr.write_arrayshape(datasetdims);
    }
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
    auto global_data = collect_global_data(h_data, xzarr_ptr->get_dimnames());
    if (my_rank == 0) {
      xzarr_ptr->write_to_array(global_data);
      xzarr_ptr->write_arrayshape(datasetdims);
    }
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
    T recv_data = data;
    if (std::same_as<T, unsigned int>) {
      MPI_Reduce(&data, &recv_data, 1, MPI_UNSIGNED, MPI_SUM, 0, comm);
    }
    if (my_rank == 0) {
      xzarr_ptr->write_to_array(recv_data);
      xzarr_ptr->write_arrayshape(datasetdims);
    }
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
    std::vector<int> distributed_sizes(comm_size);
    std::vector<int> receive_displacements(comm_size, 0), receive_counts(comm_size, 0);
    int local_size = h_data.extent(0);

    // Since there is no defined dimensions for ragged arrays collect the array sizes directly
    MPI_Gather(&local_size, 1, MPI_INT, distributed_sizes.data(), 1, MPI_INT, 0, comm);
    int global_size = std::accumulate(distributed_sizes.begin(), distributed_sizes.end(), 0);
    Kokkos::View<T *, HostSpace> global_data("global_output_data", global_size);

    // Calculate the receive counts and displacements (meaningful only for process 0)
    for (int i = 0; i < comm_size; i++) {
      receive_counts[i] = distributed_sizes[i];
      if (i > 0) receive_displacements[i] = receive_displacements[i - 1] + distributed_sizes[i - 1];
    }

    // Collect the global data
    collect_global_array(global_data.data(), h_data.data(), h_data.extent(0), receive_counts.data(),
                         receive_displacements.data());

    // IMPORTANT: In order to write the data correctly this function requires a
    // global superdroplet ordering. This is constructed based on the
    // superdroplet indices, which are unique globally. This means however, that
    // the superdroplet index observer must ALWAYS be the first in the output
    // order, so that subsequent writes use the correct ordering.
    if (my_rank == 0) {
      // If this is the superdrop index fill the ordering array
      if (std::same_as<T, unsigned int>) {
        for (size_t i = 0; i < global_data.extent(0); i++)
          global_superdroplet_ordering.get()->at(global_data[i]) = i;
      }
      // Correctly orders the data based on the global superdroplet ordering
      Kokkos::View<T *, HostSpace> global_write_data("global_write_data", global_size);
      for (size_t i = 0; i < global_data.extent(0); i++)
        global_write_data[i] = global_data[global_superdroplet_ordering.get()->at(i)];

      // Writes the final array
      xzarr.write_to_array(global_write_data);
      xzarr.write_ragged_arrayshape();
    }
  }
};

#endif  // LIBS_ZARR_COLLECTIVE_DATASET_HPP_
