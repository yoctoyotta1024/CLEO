/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: variable_metadata.hpp
 * Project: zarr
 * Created Date: Wednesday 26th November 2026
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct bundling all the metadata associated with a variable written to a Zarr array / Xarray
 * dataset.
 */

#ifndef LIBS_ZARR_VARIABLE_METADATA_HPP_
#define LIBS_ZARR_VARIABLE_METADATA_HPP_

#include <string_view>

/**
 * @brief Metadata describing a variable written to a Zarr array that is part of an Xarray dataset.
 *
 * Bundles together what's needed to write a variable's data (scale_factor, to convert CLEO's
 * dimensionless data into physical units) and what's needed to document it for a human or tool
 * reading the dataset (units, long_name, standard_name, description). Replaces passing
 * name/units/scale_factor around as three separate parameters, and gives CLEO one place to define
 * each output variable's metadata (see libs/observers/cleo_variables_metadata.hpp).
 */
struct VariableMetadata {
  std::string_view name;          /**< name of variable (key of the array in the Zarr store) */
  std::string_view units;         /**< units of variable; CF convention "1" if dimensionless */
  double scale_factor;            /**< scale factor converts dimensionless variable to "units" */
  std::string_view long_name;     /**< human-readable name of the variable */
  std::string_view standard_name; /**< CF standard_name of the variable (empty if none applies) */
  std::string_view description;   /**< free-text description of what the variable represents */

  /**
   * @brief Constructs a VariableMetadata object.
   *
   * @param name Name of variable (key of the array in the Zarr store).
   * @param units Units of variable's data; CF convention "1" if dimensionless.
   * @param scale_factor Scale factor to convert CLEO's dimensionless data to "units".
   * @param long_name Human-readable name of the variable.
   * @param standard_name CF standard_name of the variable (empty if none applies).
   * @param description Free-text description of what the variable represents.
   */
  constexpr VariableMetadata(const std::string_view name, const std::string_view units,
                             const double scale_factor, const std::string_view long_name,
                             const std::string_view standard_name,
                             const std::string_view description)
      : name(name),
        units(units),
        scale_factor(scale_factor),
        long_name(long_name),
        standard_name(standard_name),
        description(description) {}
};

#endif  // LIBS_ZARR_VARIABLE_METADATA_HPP_
