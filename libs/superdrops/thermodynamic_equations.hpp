/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: thermodynamic_equations.hpp
 * Project: superdrops
 * Created Date: Wednesday 25th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 9th February 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for functions that return Left Hand Side of thermodynamic equations. Unless stated
 * otherwise, equations referenced as (eqn [X.YY]) are from "An Introduction To Clouds From The
 * Microscale to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 * */

#ifndef LIBS_SUPERDROPS_THERMODYNAMIC_EQUATIONS_HPP_
#define LIBS_SUPERDROPS_THERMODYNAMIC_EQUATIONS_HPP_

#include <cassert>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "../cleoconstants.hpp"
#include "./superdrop.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/**
 * @brief Calculate the specific heat capacity of a moist parcel of air.
 *
 * @param qvap The vapor mass mixing ratio.
 * @param qcond The condensate mass mixing ratio.
 * @return The (dimensionless) specific heat capacity.
 */
KOKKOS_INLINE_FUNCTION
double moist_specifc_heat(const double qvap, const double qcond) {
  return dlc::Cp_dry + dlc::Cp_v * (qvap) + dlc::C_l * (qcond);
}

/**
 * @brief Calculate the supersaturation ratio given the saturation pressure, ambient pressure,
 * and vapor mass mixing ratio.
 *
 * supersaturation ratio, 's_ratio', = p_vapour/psat (i.e. is equivalent to the relative humidity)
 *
 * @param press The ambient pressure.
 * @param qvap The vapor mass mixing ratio.
 * @param psat The saturation pressure.
 * @return The supersaturation ratio.
 */
KOKKOS_INLINE_FUNCTION
double supersaturation_ratio(const double press, const double qvap, const double psat) {
  return (press * qvap) / ((dlc::Mr_ratio + qvap) * psat);
}

/**
 * @brief Calculate the Raoult and Kelvin factors for the Kohler curve.
 *
 * Calculates 1) value of 'a' in Raoult factor (exp^(a/r)) to account for effect of dissolved solute
 * on radial growth of droplet and 2) Calculates value of 'b' in Kelvin factor (1-b/r^3) to account
 * for curvature on radial growth of droplet. Equations [X.YY] are from "An Introduction To Clouds
 * From The Microscale to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 *
 * @param drop The superdroplet.
 * @param temp The ambient temperature.
 * @return A Kokkos pair containing 'a' and 'b' factors in that order.
 */
KOKKOS_INLINE_FUNCTION
Kokkos::pair<double, double> kohler_factors(const Superdrop &drop, const double temp) {
  constexpr double akoh_constant = 3.3e-7 / (dlc::TEMP0 * dlc::R0);
  const auto akoh = akoh_constant / temp;  // dimensionless version of eqn [6.24]

  constexpr double bkoh_constant = 4.3e-6 * dlc::RHO0 / dlc::MR0;
  const auto msol = drop.get_msol();
  const auto ionic = drop.get_ionic();
  const auto mr_sol = drop.get_mr_sol();
  const auto bkoh = bkoh_constant * msol * ionic / mr_sol;  // dimensionless version of eqn [6.22]

  return {akoh, bkoh};  // {a, b} = {Raoult, Kelvin} Kohler factors
}

/**
 * @brief Calculate the equilibrium vapor pressure of water over liquid water,
 * i.e. the saturation pressure.
 *
 * Equation adapted from Bjorn Steven's "make_tetens" Python function from his module
 * "moist_thermodynamics.saturation_vapour_pressures" available upon request on gitlab. Original
 * paper for formula is Murray, F. W. (1967) "On the Computation of Saturation Vapor Pressure",
 * Journal of Applied Meteorology and Climatology 6, 203–204.
 *
 * Note: function starts with conversion from dimentionless to real temperature [Kelvin],
 * TEMP = temp*Temp0, and returns dimensionless pressure from real psat = PSAT/P0.
 *
 * @param temp The (dimensionless) ambient temperature.
 * @return The (dimensionless) saturation pressure.
 */
KOKKOS_FUNCTION
double saturation_pressure(const double temp);

/**
 * @brief Calculate the equilibrium vapor pressure of water over liquid water,
 * i.e. the saturation pressure.
 *
 * Equation adapted from Python module typhon.physics.thermodynamics.e_eq_water_mk with conversion
 * to real TEMP /K = temp*Temp0 and return dimensionless psat from real psat, psat = PSAT/P0.
 *
 * @param temp The (dimensionless) temperature.
 * @return The (dimensionless) saturation pressure.
 */
KOKKOS_FUNCTION
double saturation_pressure_murphy_koop(const double temp);

/**
 * @brief Calculate the sum of the heat and vapor diffusion factors for
 * condensation-diffusion growth equation.
 *
 * Calculate the sum of heat and vapor diffusion factors 'Fkl' and 'Fdl' respectively for
 * condensation-diffusion growth equation of droplet radius. Equations [X.YY] are from "An
 * Introduction To Clouds From The Microscale to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 *
 * @param press The ambient pressure.
 * @param temp The ambient temperature.
 * @param psat The saturation pressure.
 * @return The (dimensionless) diffusion factor.
 */
KOKKOS_FUNCTION
double diffusion_factor(const double press, const double temp, const double psat);

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

/**
 * @brief Calculate the equilibrium vapor pressure of water over liquid water,
 * i.e. the saturation pressure.
 *
 * Equation adapted from Bjorn Steven's "make_tetens" python function from his module
 * "moist_thermodynamics.saturation_vapour_pressures" available upon request on gitlab. Original
 * paper for formula is Murray, F. W. (1967) "On the Computation of Saturation Vapor Pressure",
 * Journal of Applied Meteorology and Climatology 6, 203–204.
 *
 * Note: function starts with conversion from dimentionless to real temperature [Kelvin],
 * TEMP = temp*Temp0, and returns dimensionless pressure from real psat = PSAT/P0.
 *
 * @param temp The (dimensionless) ambient temperature.
 * @return The (dimensionless) saturation pressure.
 */
KOKKOS_FUNCTION
double saturation_pressure(const double temp) {
  assert((temp > 0) && "psat ERROR: temperature must be larger than 0K.");

  constexpr double A = 17.4146;     // constants from Bjorn Gitlab originally from paper
  constexpr double B = 33.639;      // ditto
  constexpr double TREF = 273.16;   // Triple point temperature [K] of water
  constexpr double PREF = 611.655;  // Triple point pressure [Pa] of water

  const auto T = double{temp * dlc::TEMP0};  // real T [K]

  return (PREF * Kokkos::exp(A * (T - TREF) / (T - B))) / dlc::P0;  // dimensionless psat
}

/**
 * @brief Calculate the equilibrium vapor pressure of water over liquid water,
 * i.e. the saturation pressure.
 *
 * Equation adapted from Python module typhon.physics.thermodynamics.e_eq_water_mk with conversion
 * to real TEMP /K = temp*Temp0 and return dimensionless psat from real psat, psat = PSAT/P0.
 *
 * @param temp The (dimensionless) temperature.
 * @return The (dimensionless) saturation pressure.
 */
KOKKOS_FUNCTION
double saturation_pressure_murphy_koop(const double temp) {
  assert((temp > 0) && "psat ERROR: temperature must be larger than 0K.");

  const auto T = double{temp * dlc::TEMP0};  // real T [K]

  const auto lnpsat =
      double{54.842763  // ln(psat) [Pa]
             - 6763.22 / T - 4.21 * log(T) + 0.000367 * T +
             tanh(0.0415 * (T - 218.8)) * (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T)};

  return Kokkos::exp(lnpsat) / dlc::P0;  // dimensionless psat
}

/**
 * @brief Calculate the sum of the heat and vapor diffusion factors for
 * condensation-diffusion growth equation.
 *
 * Calculate the sum of heat and vapor diffusion factors 'Fkl' and 'Fdl' respectively for
 * condensation-diffusion growth equation of droplet radius. Equations [X.YY] are from "An
 * Introduction To Clouds From The Microscale to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 *
 * @param press The ambient pressure.
 * @param temp The ambient temperature.
 * @param psat The saturation pressure.
 * @return The (dimensionless) diffusion factor.
 */
KOKKOS_FUNCTION
double diffusion_factor(const double press, const double temp, const double psat) {
  constexpr double A = 7.11756e-5;                             // coefficient for T^2 in T*[eq.7.24]
  constexpr double B = 4.38127686e-3;                          // coefficient for T in T*[eq.7.24]
  constexpr double LATENT_RGAS_V = DC::LATENT_V / DC::RGAS_V;  // for fkl diffusion factor calc
  constexpr double D = 4.012182971e-5;                         // constants in equation [eq.7.26]

  const auto TEMP = double{temp * dlc::TEMP0};
  const auto PRESS = double{press * dlc::P0};
  const auto PSAT = double{psat * dlc::P0};

  const auto THERMK =
      double{A * Kokkos::pow(TEMP, 2.0) + TEMP * B};  // K*TEMP with K from [eq.7.24] (for fkl)
  const auto DIFFUSE_V = double{(D / PRESS * Kokkos::pow(TEMP, 1.94)) /
                                DC::RGAS_V};  // 1/R_v * D_v from [eq 7.26] (for fdl)

  const auto fkl =
      double{(LATENT_RGAS_V / TEMP - 1.0) * DC::LATENT_V / (THERMK * dlc::F0)};  // fkl eqn [7.23]
  const auto fdl = double{TEMP / (DIFFUSE_V * PSAT) / dlc::F0};                  // fdl eqn [7.25]

  return dlc::Rho_l * (fkl + fdl);  // total constant from sum of diffusion factors
}

#endif  // LIBS_SUPERDROPS_THERMODYNAMIC_EQUATIONS_HPP_
