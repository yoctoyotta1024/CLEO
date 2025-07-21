/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: zarr_metadata.hpp
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

#ifndef LIBS_ZARR_ZARR_METADATA_HPP_
#define LIBS_ZARR_ZARR_METADATA_HPP_

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

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
std::string vec_to_string(const std::vector<size_t>& vals);

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
                                   const std::string_view dtype);

/**
 * @brief Class for generating metadata required for a Zarr array.
 *
 * This class generates compulsory metadata for the Zarr array .zarray JSON file.
 * The metadata includes information such as chunk shape, data type, etc.
 *
 * @tparam T The data type of the array's elements.
 */
template <typename T>
class ZarrMetadata {
 private:
  std::string part_zarrmetadata; /**< Metadata required for zarr array excluding array's shape */

 public:
  /**
   * @brief Constructs a ZarrMetadata object.
   *
   * Constructs a ZarrMetadata object with the given chunk shape and data type.
   *
   * @param chunkshape The shape of the chunks used to store array data.
   * @param dtype The data type of the array's elements in Zarr format (e.g., "<f8" for double).
   */
  ZarrMetadata(const std::vector<size_t>& chunkshape, const std::string_view dtype)
      : part_zarrmetadata(make_part_zarrmetadata(chunkshape, dtype)) {}

  /**
   * @brief Generates metadata for the Zarr array.
   *
   * Generates compulsory metadata for the Zarr array .zarray JSON file.
   *
   * @param arrayshape The shape of the Zarr array.
   * @return A string containing the metadata for the Zarr array.
   */
  std::string operator()(const std::vector<size_t>& arrayshape) const {
    const auto metadata = std::string(
        "{\n"
        "  \"shape\": " +
        vec_to_string(arrayshape) + ",\n" + std::string(part_zarrmetadata) + "\n}");

    return metadata;
  }
};

template <>
class ZarrMetadata<uint64_t> : public ZarrMetadata<void> {
  constexpr static char dtype[] = "<u8";

 public:
  explicit ZarrMetadata(const std::vector<size_t>& chunkshape)
      : ZarrMetadata<void>(chunkshape, dtype) {}
};

template <>
class ZarrMetadata<uint32_t> : public ZarrMetadata<void> {
  constexpr static char dtype[] = "<u4";

 public:
  explicit ZarrMetadata(const std::vector<size_t>& chunkshape)
      : ZarrMetadata<void>(chunkshape, dtype) {}
};

template <>
class ZarrMetadata<double> : public ZarrMetadata<void> {
  constexpr static char dtype[] = "<f8";

 public:
  explicit ZarrMetadata(const std::vector<size_t>& chunkshape)
      : ZarrMetadata<void>(chunkshape, dtype) {}
};

template <>
class ZarrMetadata<float> : public ZarrMetadata<void> {
  constexpr static char dtype[] = "<f4";

 public:
  explicit ZarrMetadata(const std::vector<size_t>& chunkshape)
      : ZarrMetadata<void>(chunkshape, dtype) {}
};

#endif  // LIBS_ZARR_ZARR_METADATA_HPP_
