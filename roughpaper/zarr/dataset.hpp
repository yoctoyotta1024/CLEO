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
  ZarrGroup<Store>& zarr;  ///< Reference to the zarr group object.

 public:
  /**
   * @brief Constructs a Dataset with the specified store object.
   *
   * This constructor initializes a Dataset with the provided store object by initialising a
   * ZarrGroup and writing some additional metatdata for Xarray and NetCDF.
   *
   * @param store The store object associated with the Dataset.
   */
  explicit Dataset(const Store store) : zarr(store) {
    store[".zattrs"] =
        "{\n"
        "  \"creator\": \"Clara Bayley\",\n"
        "  \"title\": \"Zarr Group for Data Output from CLEO\""
        "\n}";
  }

  void create_zarr_array() {}
};

#endif  // ROUGHPAPER_ZARR_DATASET_HPP_
