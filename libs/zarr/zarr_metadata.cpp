/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: zarr_metadata.cpp
 * Project: zarr
 * Created Date: Wednesday 22nd May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * structs to generate metadata to make zarr arrays
 */

#include "./zarr_metadata.hpp"

/**
 * @brief Converts a vector of integers into a single list written as a string.
 *
 * Given vector of a type convertible to a string with values [a, b, c, ..., z], function returns
 * the string "[a, b, c, ..., z]" with elements separated by commas and enclosed in square brackets.
 * Function is useful for converting vectors representing the shape of chunks and arrays etc. into
 * a string format for metadata json files.
 *
 * @param vals The vector values of a type convertible to a string.
 * @return A string representation of the vector.
 */
std::string vec_to_string(const std::vector<size_t>& vals) {
  auto vals_str = std::string{"["};
  for (const auto& v : vals) {
    vals_str += std::to_string(v) + ", ";
  }
  vals_str.erase(vals_str.size() - 2);  // delete last ", "
  vals_str += "]";
  return vals_str;
}

/**
 * @brief Generates part of the metadata for a Zarr array .zarray json file.
 *
 * This function constructs a string containing all the compulsory metadata of a Zarr array for its
 * .zarray json file, excluding the array's shape.
 *
 * @param chunkshape The shape of individual data chunks along each dimension.
 * @param dtype The data type stored in the arrays (e.g., "<f8").
 * @return A string view containing the partial metadata for the Zarr array.
 */
std::string make_part_zarrmetadata(const std::vector<size_t>& chunkshape,
                                   const std::string_view dtype) {
  const auto chunkshape_str = vec_to_string(chunkshape);  // shape of each chunk of array
  const auto order = 'C';  // layout of bytes in each chunk of array in storage ('C' or 'F')
  const auto compressor = std::string{"null"};  // compression of data when writing to store
  const auto fill_value = std::string{"null"};  // fill value for empty datapoints in array
  const auto filters = std::string{"null"};     // codec configurations for compression
  const auto zarr_format = '2';                 // storage spec. version 2

  const auto part_zarrmetadata = std::string("  \"chunks\": " + chunkshape_str +
                                             ",\n"
                                             "  \"dtype\": \"" +
                                             std::string(dtype) +
                                             "\",\n"
                                             "  \"order\": \"" +
                                             order +
                                             "\",\n"
                                             "  \"compressor\": " +
                                             compressor +
                                             ",\n"
                                             "  \"fill_value\": " +
                                             fill_value +
                                             ",\n"
                                             "  \"filters\": " +
                                             filters +
                                             ",\n"
                                             "  \"zarr_format\": " +
                                             zarr_format);
  return part_zarrmetadata;
}
