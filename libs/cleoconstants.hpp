/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cleoconstants.hpp
 * Project: libs
 * Created Date: Monday 29th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * header file with namespaces for (physical) constants used by CLEO
 * _Note:_ All letters in CAPITALS indicates constants have dimensions
 * (or are naturally dimensionless as for IONIC). Only first letter
 * capitalised indicates value is dimensionless (version of a) constant.
 */

#ifndef LIBS_CLEOCONSTANTS_HPP_
#define LIBS_CLEOCONSTANTS_HPP_

#include <cstdint>
#include <limits>
#include <numbers>

/**
 * @brief Namespace containing values of constants with dimensions.
 */
namespace dimmed_constants {
constexpr double G = 9.80665;                   /**< Acceleration due to gravity [m/s^2]. */
constexpr double RGAS_UNIV = 8.314462618;       /**< Universal molar gas constant [J/Kg/K]. */
constexpr double MR_WATER = 0.01801528;         /**< Molecular mass of water [Kg/mol]. */
constexpr double MR_DRY = 0.028966216;          /**< Molecular mass of dry air [Kg/mol]. */
constexpr double RGAS_DRY = RGAS_UNIV / MR_DRY; /**< Specific gas constant for dry air [J/Kg/K]. */
constexpr double RGAS_V = RGAS_UNIV / MR_WATER; /**< Specific gas constant for water [J/Kg/K]. */

constexpr double LATENT_V = 2500930;
/**< Specific latent heat of vapourisation of water [J/Kg]  (IAPWS97 at 273.15K). */
constexpr double CP_DRY = 1004.64;
/**< Specific heat capacity (dry) air at constant pressure [J/Kg/K] ~1.400*cv_dry (ICON). */
constexpr double CP_V = 1865.01;
/**< Specific heat capacity of water vapour [J/Kg/K] (IAPWS97 at 273.15K). */
constexpr double C_L = 4192.664;
/**< Specific heat capacity of liquid water[J/Kg/K] (ICON c_l = (3.1733 + 1.0) * cp_dry). */

constexpr double RHO_DRY = 1.177; /**< Density of dry air [Kg/m^3] (at 300K). */
constexpr double RHO_L = 998.203;
/**< Density of liquid water condensing [kg/m^3] (water at 293K from SCALE-SDM). */
constexpr double DYNVISC = 18.45 * 1e-6; /**< dynamic viscosity of air [Pa s] (at 300K). */

constexpr double RHO_SOL = 2016.5; /**< Density of (dry) areosol [Kg/m^3] (NaCl from SCALE-SDM). */
constexpr double MR_SOL = 0.05844277;
/**< Molecular mass of areosol [Kg/mol] (NaCl=0.058 from SCALE-SDM). */
constexpr int IONIC = 2; /**< Degree ionic dissociation (van't Hoff factor) [dimensionless]. */

constexpr double SURFSIGMA = 7.28e-2; /**< Surface tension of water [J/m^-2]. */
}  // namespace dimmed_constants

/**
 * @brief Namespace containing constants for using characterstic scales of time, velocity,
 * temperature, ... (TIME0, TEMP0, P0, ...) etc. in order to make variables dimensionless.
 * Namespace also includes the dimensionless equivalents of some members of the
 * dimmed_constants namespace.
 *
 */
namespace dimless_constants {

namespace DC = dimmed_constants;

/* characterstic scales */
constexpr double W0 = 1.0;                        /**< Characteristic velocity [m/s]. */
constexpr double TIME0 = 1000.0;                  /**< Timescale [s]. */
constexpr double COORD0 = TIME0 * W0;             /**< Coordinate grid scale [m]. */
constexpr double VOL0 = COORD0 * COORD0 * COORD0; /**< Volume scale [m^3]. */

constexpr double CP0 = DC::CP_DRY; /**< Heat capacity [J/Kg/K]. */
constexpr double MR0 = DC::MR_DRY; /**< Molecular molar mass [Kg/mol]. */
constexpr double R0 = 1e-6;        /**< Droplet radius lengthscale [m]. */

constexpr double P0 = 100000.0;                 /**< Pressure [Pa]. */
constexpr double TEMP0 = 273.15;                /**< Temperature [K]. */
constexpr double RHO0 = P0 / (CP0 * TEMP0);     /**< Density [Kg/m^3]. */
constexpr double MASS0 = R0 * R0 * R0 * RHO0;   /**< Mass [Kg]. */
constexpr double MASS0grams = MASS0 * 1000;     /**< Mass [g]. */
constexpr double F0 = TIME0 / (RHO0 * R0 * R0); /**< Droplet condensation-diffusion factors. */

/* dimensionless constants */
constexpr double Mr_ratio = DC::MR_WATER / DC::MR_DRY; /**< Dimensionless molecular mass ratio. */
constexpr double Cp_dry = DC::CP_DRY / CP0;            /**< Dimensionless CP_DRY. */
constexpr double Cp_v = DC::CP_V / CP0;                /**< Dimensionless CP_V. */
constexpr double C_l = DC::C_L / CP0;                  /**< Dimensionless C_L. */
constexpr double Latent_v = DC::LATENT_V / (TEMP0 * CP0); /**< Dimensionless LATENT_V. */
constexpr double Rgas_dry = DC::RGAS_DRY / CP0;           /**< Dimensionless RGAS_DRY. */
constexpr double Rgas_v = DC::RGAS_V / CP0;               /**< Dimensionless RGAS_V. */
constexpr double Rho_dry = DC::RHO_DRY / RHO0;            /**< Dimensionless RHO_DRY. */
constexpr double Rho_l = DC::RHO_L / RHO0;                /**< Dimensionless RHO_L. */
constexpr double Rho_sol = DC::RHO_SOL / RHO0;            /**< Dimensionless RHO_SOL. */
constexpr double Mr_sol = DC::MR_SOL / MR0;               /**< Dimensionless MR_SOL. */
constexpr int IONIC = DC::IONIC;                          /**< Dimensionless IONIC. */

constexpr double surfconst = 4.0 * DC::SURFSIGMA * std::numbers::pi * R0 * R0;
/**< Surface tension energy constant (see collisionkinetics) */
}  // namespace dimless_constants

/**
 * @brief max/min values for certain types, e.g. to use values from C++ standard numeric limits on
 * GPUs.
 *
 */
namespace LIMITVALUES {
constexpr unsigned int uintmax = std::numeric_limits<unsigned int>::max(); /**< Max unsigned int. */
constexpr uint64_t uint64_t_max = std::numeric_limits<uint64_t>::max(); /**< Max 64 byte u-int. */

constexpr unsigned int oob_gbxindex = uintmax; /**< value for an out of bounds gbxindex */
constexpr double llim = -1.0 * std::numeric_limits<double>::max(); /**< Maximum negative double. */
constexpr double ulim = std::numeric_limits<double>::max();        /**< Maximum positive double. */
}  // namespace LIMITVALUES

#endif  // LIBS_CLEOCONSTANTS_HPP_
