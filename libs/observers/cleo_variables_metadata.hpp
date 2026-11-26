/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cleo_variables_metadata.hpp
 * Project: observers
 * Created Date: Wednesday 26th November 2026
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Centralised definitions of the metadata (name, units, scale_factor, long_name, standard_name,
 * description) for every variable CLEO's observers can write to a Zarr dataset.
 */

#ifndef LIBS_OBSERVERS_CLEO_VARIABLES_METADATA_HPP_
#define LIBS_OBSERVERS_CLEO_VARIABLES_METADATA_HPP_

#include "../cleoconstants.hpp"
#include "zarr/variable_metadata.hpp"

namespace dlc = dimless_constants;

/**
 * @brief Namespace holding the single, centralised definition of the metadata (name, units,
 * scale_factor, long_name, standard_name, description) for every variable CLEO's observers can
 * write to a Zarr dataset.
 *
 * Each observer that writes a variable should use the corresponding VariableMetadata constant
 * defined here instead of hardcoding its own copy of the variable's name/units/scale_factor, so
 * there is exactly one place to look up, or correct, a variable's metadata.
 */
namespace cleo_metadata {

inline constexpr VariableMetadata time{"time", "s",    dlc::TIME0,
                                       "time", "time", "time elapsed since start of simulation"};

inline constexpr VariableMetadata gbxindex{
    "gbxindex",      "1", 1.0,
    "gridbox index", "",  "index of each gridbox in the domain (unique gridbox ID)"};

inline constexpr VariableMetadata press{
    "press", "hPa", dlc::P0 / 100, "air pressure", "air_pressure", "pressure at centre of gridbox"};

inline constexpr VariableMetadata temp{"temp",
                                       "K",
                                       dlc::TEMP0,
                                       "air temperature",
                                       "air_temperature",
                                       "temperature at centre of gridbox"};

inline constexpr VariableMetadata qvap{"qvap",
                                       "g/kg",
                                       1000.0,
                                       "vapour mass mixing ratio",
                                       "humidity_mixing_ratio",
                                       "mass of water vapour per unit mass of dry air in gridbox"};

inline constexpr VariableMetadata qcond{
    "qcond", "g/kg",
    1000.0,  "liquid mass mixing ratio",
    "",      "mass of condensed (total liquid) water per unit mass of dry air in gridbox"};

inline constexpr VariableMetadata wvel{
    "wvel",
    "m/s",
    dlc::W0,
    "vertical wind velocity",
    "upward_air_velocity",
    "wind velocity in the vertical (z) direction at centre of gridbox"};

inline constexpr VariableMetadata uvel{
    "uvel",          "m/s",
    dlc::W0,         "eastward wind velocity",
    "eastward_wind", "wind velocity in the eastward (x) direction at centre of gridbox"};

inline constexpr VariableMetadata vvel{
    "vvel",           "m/s",
    dlc::W0,          "northward wind velocity",
    "northward_wind", "wind velocity in the northward (y) direction at centre of gridbox"};

inline constexpr VariableMetadata nsupers{
    "nsupers", "1", 1.0, "number of superdroplets", "", "number of superdroplets in each gridbox"};

inline constexpr VariableMetadata totnsupers{
    "totnsupers", "1",
    1.0,          "total number of superdroplets",
    "",           "total number of superdroplets in the entire domain"};

inline constexpr VariableMetadata massmom0{
    "massmom0",
    "1",
    1.0,
    "0th mass moment",
    "",
    "0th mass moment of the droplet mass distribution in each gridbox "
    "(i.e. the number of droplets)"};
inline constexpr VariableMetadata massmom0_raindrops{
    "massmom0_raindrops",
    "1",
    1.0,
    "0th mass moment of raindrops",
    "",
    "0th mass moment of the raindrop (>= 40micron) mass distribution in each gridbox "
    "(i.e. the number of droplets >= 40micron)"};

inline constexpr VariableMetadata massmom1{
    "massmom1",
    "g",
    dlc::MASS0grams,
    "1st mass moment",
    "",
    "1st mass moment of the droplet mass distribution in each gridbox "
    "(i.e. the total droplet mass)"};
inline constexpr VariableMetadata massmom1_raindrops{
    "massmom1_raindrops",
    "g",
    dlc::MASS0grams,
    "1st mass moment of raindrops",
    "",
    "1st mass moment of the raindrop (>= 40micron) mass distribution "
    "in each gridbox (i.e. total droplet mass of droplets >= 40micron)"};

inline constexpr VariableMetadata massmom2{
    "massmom2",
    "g^2",
    dlc::MASS0grams* dlc::MASS0grams,
    "2nd mass moment",
    "",
    "2nd mass moment of the droplet mass distribution in each gridbox "
    "(proportional to radar reflectivity in Rayleigh regime)"};
inline constexpr VariableMetadata massmom2_raindrops{
    "massmom2_raindrops",
    "g^2",
    dlc::MASS0grams* dlc::MASS0grams,
    "2nd mass moment of raindrops",
    "",
    "2nd mass moment of the raindrop (>= 40micron) mass distribution in each gridbox "
    "(proportional to radar reflectivity in Rayleigh regime but "
    "accounting for only droplets >= 40micron)"};

inline constexpr VariableMetadata delta_massmom0_microphys{
    "delta_massmom0_microphys",
    "1",
    1.0,
    "change in 0th mass moment from microphysics",
    "",
    "change in the 0th mass moment of the droplet distribution due to microphysics"};
inline constexpr VariableMetadata delta_massmom1_microphys{
    "delta_massmom1_microphys",
    "g",
    dlc::MASS0grams,
    "change in 1st mass moment from microphysics",
    "",
    "change in the 1st mass moment of the droplet distribution due to microphysics"};
inline constexpr VariableMetadata delta_massmom2_microphys{
    "delta_massmom2_microphys",
    "g^2",
    dlc::MASS0grams* dlc::MASS0grams,
    "change in 2nd mass moment from microphysics",
    "",
    "change in the 2nd mass moment of the droplet distribution due to microphysics"};

inline constexpr VariableMetadata delta_massmom0_motion{
    "delta_massmom0_motion",
    "1",
    1.0,
    "change in 0th mass moment from motion",
    "",
    "change in the 0th mass moment of the droplet distribution due to superdroplet motion"};
inline constexpr VariableMetadata delta_massmom1_motion{
    "delta_massmom1_motion",
    "g",
    dlc::MASS0grams,
    "change in 1st mass moment from motion",
    "",
    "change in the 1st mass moment of the droplet distribution due to superdroplet motion"};
inline constexpr VariableMetadata delta_massmom2_motion{
    "delta_massmom2_motion",
    "g^2",
    dlc::MASS0grams* dlc::MASS0grams,
    "change in 2nd mass moment from motion",
    "",
    "change in the 2nd mass moment of the droplet distribution due to superdroplet motion"};

inline constexpr VariableMetadata delta_massmom0_raindrops_microphys{
    "delta_massmom0_raindrops_microphys",
    "1",
    1.0,
    "change in 0th mass moment of raindrops from microphysics",
    "",
    "change in the 0th mass moment of the raindrop (>= 40micron) distribution "
    "due to microphysics"};
inline constexpr VariableMetadata delta_massmom1_raindrops_microphys{
    "delta_massmom1_raindrops_microphys",
    "g",
    dlc::MASS0grams,
    "change in 1st mass moment of raindrops from microphysics",
    "",
    "change in the 1st mass moment of the raindrop (>= 40micron) distribution "
    "due to microphysics"};
inline constexpr VariableMetadata delta_massmom2_raindrops_microphys{
    "delta_massmom2_raindrops_microphys",
    "g^2",
    dlc::MASS0grams* dlc::MASS0grams,
    "change in 2nd mass moment of raindrops from microphysics",
    "",
    "change in the 2nd mass moment of the raindrop (>= 40micron) distribution "
    "due to microphysics in each gridbox"};

inline constexpr VariableMetadata delta_massmom0_raindrops_motion{
    "delta_massmom0_raindrops_motion",
    "1",
    1.0,
    "change in 0th mass moment of raindrops from motion",
    "",
    "change in the 0th mass moment of the raindrop (>= 40micron) distribution "
    "due to superdroplet motion"};
inline constexpr VariableMetadata delta_massmom1_raindrops_motion{
    "delta_massmom1_raindrops_motion",
    "g",
    dlc::MASS0grams,
    "change in 1st mass moment of raindrops from motion",
    "",
    "change in the 1st mass moment of the raindrop (>= 40micron) distribution "
    "due to superdroplet motion in each gridbox"};
inline constexpr VariableMetadata delta_massmom2_raindrops_motion{
    "delta_massmom2_raindrops_motion",
    "g^2",
    dlc::MASS0grams* dlc::MASS0grams,
    "change in 2nd mass moment of raindrops from motion",
    "",
    "change in the 2nd mass moment of the raindrop (>= 40micron) distribution "
    "due to superdroplet motion in each gridbox"};

inline constexpr VariableMetadata massdelta_cond{
    "massdelta_cond",
    "g",
    dlc::MASS0grams,
    "condensation/evaporation mass change",
    "",
    "change in mass of water due to condensation/evaporation in each gridbox"};

inline constexpr VariableMetadata precip{
    "precip",
    "m",
    dlc::R0* dlc::R0* dlc::R0 / dlc::COORD0 / dlc::COORD0,
    "precipitation",
    "",
    "accumulated precipitation (equivalent to water depth) from droplets "
    "crossing the lower boundary of each gridbox"};

inline constexpr VariableMetadata sdgbxindex{
    "sdgbxindex", "1",
    1.0,          "superdroplet gridbox index",
    "",           "index of the gridbox each superdroplet occupies"};

inline constexpr VariableMetadata sdId{
    "sdId", "1", 1.0, "superdroplet ID", "", "unique identity of each superdroplet"};

inline constexpr VariableMetadata xi{
    "xi", "1", 1.0, "multiplicity", "", "number of 'real' droplets a superdroplet represents"};

inline constexpr VariableMetadata radius{
    "radius",
    "um",
    dlc::R0 * 1e6,
    "radius",
    "",
    "size of droplet assuming its spherical and accounting for volume of its (dry) solute"};

inline constexpr VariableMetadata msol{
    "msol",           "g", dlc::MASS0grams,
    "mass of solute", "",  "mass of (dry) solute in each superdroplet"};

inline constexpr VariableMetadata coord1{
    "coord1", "m", dlc::COORD0, "eastwards (x) coordinate", "", "x position of each superdroplet"};
inline constexpr VariableMetadata coord2{
    "coord2", "m", dlc::COORD0, "northwards (y) coordinate", "", "y position of each superdroplet"};
inline constexpr VariableMetadata coord3{
    "coord3", "m", dlc::COORD0, "vertical (z) coordinate", "", "z position of each superdroplet"};

inline constexpr VariableMetadata raggedcount{
    "raggedcount",
    "1",
    1.0,
    "ragged count",
    "",
    "number of superdroplets in each gridbox at each timestep (count variable for ragged arrays)"};

}  // namespace cleo_metadata

#endif  // LIBS_OBSERVERS_CLEO_VARIABLES_METADATA_HPP_
