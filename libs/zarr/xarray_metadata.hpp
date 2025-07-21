/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: xarray_metadata.hpp
 * Project: zarr
 * Created Date: Wednesday 22nd May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * structs to generate metadata to make zarr arrays conform with xarray dataset
 */

#ifndef LIBS_ZARR_XARRAY_METADATA_HPP_
#define LIBS_ZARR_XARRAY_METADATA_HPP_

#include <cassert>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make zarr array
 * compatible with Xarray and NetCDF for a certain data type, T.
 *
 * @tparam T The data type of the array.
 * @param units The units of the array's coordinates.
 * @param scale_factor The scale factor of data.
 * @param dimnames The names of each dimension of the array.
 * @return A string representing the metadata.
 */
template <typename T>
std::string xarray_metadata(const std::string_view units, const double scale_factor,
                            const std::vector<std::string>& dimnames);

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make array for
 * a raggedcount variable in Zarr compatible with Xarray and NetCDF for a certain data type, T.
 *
 * @tparam T The data type of the array.
 * @param units The units of the array's coordinates.
 * @param scale_factor The scale factor of data.
 * @param dimnames The names of each dimension of the array.
 * @param sampledimname The name of the dimension the ragged count samples.
 * @return A string representing the metadata.
 */
template <typename T>
std::string xarray_metadata(const std::string_view units, const double scale_factor,
                            const std::vector<std::string>& dimnames,
                            const std::string_view sampledimname);

template <>
std::string xarray_metadata<double>(const std::string_view units, const double scale_factor,
                                    const std::vector<std::string>& dimnames);

template <>
std::string xarray_metadata<float>(const std::string_view units, const double scale_factor,
                                   const std::vector<std::string>& dimnames);

template <>
std::string xarray_metadata<uint32_t>(const std::string_view units, const double scale_factor,
                                      const std::vector<std::string>& dimnames);

template <>
std::string xarray_metadata<uint64_t>(const std::string_view units, const double scale_factor,
                                      const std::vector<std::string>& dimnames);

template <>
std::string xarray_metadata<double>(const std::string_view units, const double scale_factor,
                                    const std::vector<std::string>& dimnames,
                                    const std::string_view sampledimname);

template <>
std::string xarray_metadata<float>(const std::string_view units, const double scale_factor,
                                   const std::vector<std::string>& dimnames,
                                   const std::string_view sampledimname);

template <>
std::string xarray_metadata<uint32_t>(const std::string_view units, const double scale_factor,
                                      const std::vector<std::string>& dimnames,
                                      const std::string_view sampledimname);
template <>
std::string xarray_metadata<uint64_t>(const std::string_view units, const double scale_factor,
                                      const std::vector<std::string>& dimnames,
                                      const std::string_view sampledimname);

#endif  // LIBS_ZARR_XARRAY_METADATA_HPP_
