/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: xarray_metadata.cpp
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

#include "./xarray_metadata.hpp"

/**
 * @brief Converts vector of strings, e.g. for names of dimensions,into a single list written
 * as a string.
 *
 * @param dims The vector of strings to be converted.
 * @return A string representing the converted list.
 */
std::string vecstr_to_string(const std::vector<std::string>& dims) {
  auto dims_str = std::string{"["};
  for (const auto& d : dims) {
    dims_str += "\"" + d + "\",";
  }
  dims_str.pop_back();  // delete last ","
  dims_str += "]";
  return dims_str;
}

/**
 * @brief Make string of with set precision out of of scale factor double.
 *
 * Use precision of limit for float (~6-7 decimal places).
 *
 * @param scale_factor The scale factor of data.
 * @return A string representing the scale factor.
 */
std::string scale_factor_string(const double scale_factor) {
  const int prec = std::numeric_limits<float>::digits10;  // precision (no. decimal digits) of float
  std::ostringstream oss;
  oss << std::scientific << std::setprecision(prec) << scale_factor;
  std::string scale_factor_str = oss.str();

  return scale_factor_str;
}

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make zarr array
 * compatible with Xarray and NetCDF.
 *
 * Metadata includes scale_factor so is only valid for floating point types (dtype "f") or
 * complex floating point (dtype "c").
 *
 * @param units The units of the array's coordinates.
 * @param scale_factor The scale factor of data.
 * @param dimnames The names of each dimension of the array.
 * @return A string representing the metadata.
 */
std::string xarray_metadata_for_floats(const std::string_view units, const double scale_factor,
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
      scale_factor_string(scale_factor) +  // scale_factor of data
      "\n}");
  return zattrs;
}

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make zarr array
 * compatible with Xarray and NetCDF.
 *
 * Metadata for integer types (e.g. dtype "u") doesn't include scale_factor so assertion
 * checks scale_factor is equal to 1.0.
 *
 * @param units The units of the array's coordinates.
 * @param scale_factor The scale factor of data, must equal 1.0.
 * @param dimnames The names of each dimension of the array.
 * @return A string representing the metadata.
 */
std::string xarray_metadata_for_ints(const std::string_view units, const double scale_factor,
                                     const std::vector<std::string>& dimnames) {
  assert((scale_factor == 1.0) && "scale_factor cannot be used on non-floating point type");
  const auto zattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dimnames) +  // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(units) + "\"" +  // units of coordinate being stored
      "\n}");
  return zattrs;
}

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make array for
 * a raggedcount variable in Zarr compatible with Xarray and NetCDF.
 *
 * Metadata includes scale_factor so is only valid for floating point types (dtype "f") or
 * complex floating point (dtype "c").
 *
 * @param units The units of the array's coordinates.
 * @param scale_factor The scale factor of data.
 * @param dimnames The names of each dimension of the array.
 * @param sampledimname The name of the dimension the ragged count samples.
 * @return A string representing the metadata.
 */
std::string raggedarray_xarray_metadata_for_floats(const std::string_view units,
                                                   const double scale_factor,
                                                   const std::vector<std::string>& dimnames,
                                                   const std::string_view sampledimname) {
  const auto zattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dimnames) +  // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(units) + "\"" +  // units of coordinate being stored
      ",\n"
      "  \"scale_factor\": " +
      scale_factor_string(scale_factor) +  // scale_factor of data
      ",\n"
      "  \"sample_dimension\": " +
      "\"" + std::string(sampledimname) + "\"" +  // name of sample dimension
      "\n}");
  return zattrs;
}

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make array for
 * a raggedcount variable in Zarr compatible with Xarray and NetCDF.
 *
 * Metadata for integer types (e.g. dtype "u") doesn't include scale_factor so assertion
 * checks scale_factor is equal to 1.0.
 *
 * @param units The units of the array's coordinates.
 * @param scale_factor The scale factor of data.
 * @param dimnames The names of each dimension of the array.
 * @param sampledimname The name of the dimension the ragged count samples.
 * @return A string representing the metadata.
 */
std::string raggedarray_xarray_metadata_for_ints(const std::string_view units,
                                                 const double scale_factor,
                                                 const std::vector<std::string>& dimnames,
                                                 const std::string_view sampledimname) {
  assert((scale_factor == 1.0) && "scale_factor cannot be used on non-floating point type");
  const auto zattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dimnames) +  // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(units) + "\"" +  // units of coordinate being stored
      ",\n"
      "  \"sample_dimension\": " +
      "\"" + std::string(sampledimname) + "\"" +  // name of sample dimension
      "\n}");
  return zattrs;
}

template <>
std::string xarray_metadata<double>(const std::string_view units, const double scale_factor,
                                    const std::vector<std::string>& dimnames) {
  return xarray_metadata_for_floats(units, scale_factor, dimnames);
}

template <>
std::string xarray_metadata<float>(const std::string_view units, const double scale_factor,
                                   const std::vector<std::string>& dimnames) {
  return xarray_metadata_for_floats(units, scale_factor, dimnames);
}

template <>
std::string xarray_metadata<uint32_t>(const std::string_view units, const double scale_factor,
                                      const std::vector<std::string>& dimnames) {
  return xarray_metadata_for_ints(units, scale_factor, dimnames);
}

template <>
std::string xarray_metadata<uint64_t>(const std::string_view units, const double scale_factor,
                                      const std::vector<std::string>& dimnames) {
  return xarray_metadata_for_ints(units, scale_factor, dimnames);
}

template <>
std::string xarray_metadata<double>(const std::string_view units, const double scale_factor,
                                    const std::vector<std::string>& dimnames,
                                    const std::string_view sampledimname) {
  return raggedarray_xarray_metadata_for_floats(units, scale_factor, dimnames, sampledimname);
}

template <>
std::string xarray_metadata<float>(const std::string_view units, const double scale_factor,
                                   const std::vector<std::string>& dimnames,
                                   const std::string_view sampledimname) {
  return raggedarray_xarray_metadata_for_floats(units, scale_factor, dimnames, sampledimname);
}

template <>
std::string xarray_metadata<uint32_t>(const std::string_view units, const double scale_factor,
                                      const std::vector<std::string>& dimnames,
                                      const std::string_view sampledimname) {
  return raggedarray_xarray_metadata_for_ints(units, scale_factor, dimnames, sampledimname);
}

template <>
std::string xarray_metadata<uint64_t>(const std::string_view units, const double scale_factor,
                                      const std::vector<std::string>& dimnames,
                                      const std::string_view sampledimname) {
  return raggedarray_xarray_metadata_for_ints(units, scale_factor, dimnames, sampledimname);
}
