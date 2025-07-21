/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: zarr_group.hpp
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

#ifndef LIBS_ZARR_ZARR_GROUP_HPP_
#define LIBS_ZARR_ZARR_GROUP_HPP_

#include <Kokkos_Core.hpp>
#include <string>

/**
 * @brief A class representing a Zarr group (i.e. collection of Zarr arrays) in a storage system.
 *
 * This class provides functionality to create a group of arrays obeying the Zarr storage
 * specification version 2 (https://zarr.readthedocs.io/en/stable/spec/v2.html) within a store
 * object that manages the storage and retrieval of data and metadata.
 *
 * @tparam Store The type of the store object used by the Zarr group.
 */
template <typename Store>
struct ZarrGroup {
 public:
  Store& store; /**< Reference to the store object. */

  /**
   * @brief Constructs a ZarrGroup with the specified store object.
   *
   * This constructor initializes a ZarrGroup with the provided store object.
   * It also writes the compulsory metatdata for the group in order to obey the Zarr storage
   * specification version 2 (https://zarr.readthedocs.io/en/stable/spec/v2.html).
   *
   * @param store The store object associated with the Zarr group.
   */
  explicit ZarrGroup(Store& store) : store(store) {
    const std::string zarr_format("2");  // storage specification version 2
    const std::string zgroupjson("{\n  \"zarr_format\": " + zarr_format + "\n}");
    store[".zgroup"] = zgroupjson;
  }
};

#endif  // LIBS_ZARR_ZARR_GROUP_HPP_
