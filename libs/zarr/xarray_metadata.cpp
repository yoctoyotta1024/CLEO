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
 * @brief Returns a ",\n  \"key\": \"value\"" .zattrs entry if value is non-empty, otherwise an
 * empty string.
 *
 * Used so optional metadata fields (long_name, standard_name, description) are only written to
 * .zattrs when they have actually been set for a variable, rather than cluttering every array's
 * attributes with blank placeholders.
 *
 * @param key The name of the .zattrs key to write (e.g. "long_name").
 * @param value The value of the key; if empty, nothing is appended.
 * @return A string with the extra .zattrs entry, or an empty string if value is empty.
 */
std::string optional_attr(const std::string_view key, const std::string_view value) {
  if (value.empty()) {
    return std::string{};
  }
  return ",\n  \"" + std::string(key) + "\": \"" + std::string(value) + "\"";
}

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make zarr array
 * compatible with Xarray and NetCDF.
 *
 * Metadata includes scale_factor so is only valid for floating point types (dtype "f") or
 * complex floating point (dtype "c").
 *
 * @param meta The variable's metadata (units, scale_factor, long_name, standard_name, description).
 * @param dimnames The names of each dimension of the array.
 * @return A string representing the metadata.
 */
std::string xarray_metadata_for_floats(const VariableMetadata& meta,
                                       const std::vector<std::string>& dimnames) {
  const auto zattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dimnames) +  // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(meta.units) + "\"" +  // units of coordinate being stored
      ",\n"
      "  \"scale_factor\": " +
      scale_factor_string(meta.scale_factor) +  // scale_factor of data
      optional_attr("long_name", meta.long_name) +
      optional_attr("standard_name", meta.standard_name) +
      optional_attr("description", meta.description) +
      "\n}");
  return zattrs;
}

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make zarr array
 * compatible with Xarray and NetCDF.
 *
 * Metadata for integer types (e.g. dtype "u") doesn't include scale_factor so error-test included
 * to check scale_factor is equal to 1.0.
 *
 * @param meta The variable's metadata (units, scale_factor, long_name, standard_name, description).
 * @param dimnames The names of each dimension of the array.
 * @return A string representing the metadata.
 */
std::string xarray_metadata_for_ints(const VariableMetadata& meta,
                                     const std::vector<std::string>& dimnames) {
  if (meta.scale_factor != 1.0) {
    throw std::runtime_error("scale_factor cannot be used on non-floating point type");
  }
  const auto zattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dimnames) +  // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(meta.units) + "\"" +  // units of coordinate being stored
      optional_attr("long_name", meta.long_name) +
      optional_attr("standard_name", meta.standard_name) +
      optional_attr("description", meta.description) +
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
 * @param meta The variable's metadata (units, scale_factor, long_name, standard_name, description).
 * @param dimnames The names of each dimension of the array.
 * @param sampledimname The name of the dimension the ragged count samples.
 * @return A string representing the metadata.
 */
std::string raggedarray_xarray_metadata_for_floats(const VariableMetadata& meta,
                                                   const std::vector<std::string>& dimnames,
                                                   const std::string_view sampledimname) {
  const auto zattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dimnames) +  // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(meta.units) + "\"" +  // units of coordinate being stored
      ",\n"
      "  \"scale_factor\": " +
      scale_factor_string(meta.scale_factor) +  // scale_factor of data
      ",\n"
      "  \"sample_dimension\": " +
      "\"" + std::string(sampledimname) + "\"" +  // name of sample dimension
      optional_attr("long_name", meta.long_name) +
      optional_attr("standard_name", meta.standard_name) +
      optional_attr("description", meta.description) +
      "\n}");
  return zattrs;
}

/**
 * @brief Make string of array attributes metadata for .zattrs json which is used to make array for
 * a raggedcount variable in Zarr compatible with Xarray and NetCDF.
 *
 * Metadata for integer types (e.g. dtype "u") doesn't include scale_factor so error-test included
 * to check scale_factor is equal to 1.0.
 *
 * @param meta The variable's metadata (units, scale_factor, long_name, standard_name, description).
 * @param dimnames The names of each dimension of the array.
 * @param sampledimname The name of the dimension the ragged count samples.
 * @return A string representing the metadata.
 */
std::string raggedarray_xarray_metadata_for_ints(const VariableMetadata& meta,
                                                 const std::vector<std::string>& dimnames,
                                                 const std::string_view sampledimname) {
  if (meta.scale_factor != 1.0) {
    throw std::runtime_error("scale_factor cannot be used on non-floating point type");
  }
  const auto zattrs = std::string(
      "{\n"
      "  \"_ARRAY_DIMENSIONS\": " +
      vecstr_to_string(dimnames) +  // names of each dimension of array
      ",\n"
      "  \"units\": " +
      "\"" + std::string(meta.units) + "\"" +  // units of coordinate being stored
      ",\n"
      "  \"sample_dimension\": " +
      "\"" + std::string(sampledimname) + "\"" +  // name of sample dimension
      optional_attr("long_name", meta.long_name) +
      optional_attr("standard_name", meta.standard_name) +
      optional_attr("description", meta.description) +
      "\n}");
  return zattrs;
}

template <>
std::string xarray_metadata<double>(const VariableMetadata& meta,
                                    const std::vector<std::string>& dimnames) {
  return xarray_metadata_for_floats(meta, dimnames);
}

template <>
std::string xarray_metadata<float>(const VariableMetadata& meta,
                                   const std::vector<std::string>& dimnames) {
  return xarray_metadata_for_floats(meta, dimnames);
}

template <>
std::string xarray_metadata<uint32_t>(const VariableMetadata& meta,
                                      const std::vector<std::string>& dimnames) {
  return xarray_metadata_for_ints(meta, dimnames);
}

template <>
std::string xarray_metadata<uint64_t>(const VariableMetadata& meta,
                                      const std::vector<std::string>& dimnames) {
  return xarray_metadata_for_ints(meta, dimnames);
}

template <>
std::string xarray_metadata<double>(const VariableMetadata& meta,
                                    const std::vector<std::string>& dimnames,
                                    const std::string_view sampledimname) {
  return raggedarray_xarray_metadata_for_floats(meta, dimnames, sampledimname);
}

template <>
std::string xarray_metadata<float>(const VariableMetadata& meta,
                                   const std::vector<std::string>& dimnames,
                                   const std::string_view sampledimname) {
  return raggedarray_xarray_metadata_for_floats(meta, dimnames, sampledimname);
}

template <>
std::string xarray_metadata<uint32_t>(const VariableMetadata& meta,
                                      const std::vector<std::string>& dimnames,
                                      const std::string_view sampledimname) {
  return raggedarray_xarray_metadata_for_ints(meta, dimnames, sampledimname);
}

template <>
std::string xarray_metadata<uint64_t>(const VariableMetadata& meta,
                                      const std::vector<std::string>& dimnames,
                                      const std::string_view sampledimname) {
  return raggedarray_xarray_metadata_for_ints(meta, dimnames, sampledimname);
}
