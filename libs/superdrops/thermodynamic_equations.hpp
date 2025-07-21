/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: thermodynamic_equations.hpp
 * Project: superdrops
 * Created Date: Wednesday 25th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors: Florian Poydenot
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for functions that return Left Hand Side of thermodynamic equations. Unless stated
 * otherwise, equations referenced as (eqn [X.YY]) are from "An Introduction To Clouds From The
 * Microscale to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 */

#ifndef LIBS_SUPERDROPS_THERMODYNAMIC_EQUATIONS_HPP_
#define LIBS_SUPERDROPS_THERMODYNAMIC_EQUATIONS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <cassert>

#include "../cleoconstants.hpp"
#include "superdrop.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/**
 * @brief Calculate the specific heat capacity of moist air.
 *
 * This function calculates the specific heat capacity of a moist parcel of air using the specific
 * heat of dry air, the specific heat of water vapor, and the specific heat of condensed water, and
 * the vapour and liquid mass mixing ratios for that parcel of air.
 *
 * @param qvap The vapor mass mixing ratio.
 * @param qcond The liquid (condensate) mass mixing ratio.
 * @return The specific heat capacity of the air parcel.
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
 * _Note:_ Function starts with conversion from dimentionless to real temperature [Kelvin],
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

/**
 * @brief Calculate the ventilation factor for the condensation-diffusion growth equation.
 *
 * Equation for ventilation factor, $f_v$, is fit to data from Kinzer and Gunn (1951) and from
 * Pruppacher and Rasmussen (1979) according to Florian Poydenot, whereby
 * \f$ f_v = 1 + \frac{1}{\frac{1}{c_1R^\alpha} + \frac{1}{c_2R^\beta}} \f$
 * where \f$ c_1 = 6.954*10^7 \f$, \f$ \alpha=1.963 \f$, \f$ c_2=1.069*10^3 \f$,
 * \f$ \beta=0.702 \f$, and \f$R\f$ is the radius
 * of the water droplet in [m].
 *
 * @param radius The droplet radius.
 * @return The (dimensionless) ventilation factor.
 */
KOKKOS_FUNCTION
double ventilation_factor(const double radius);

#endif  // LIBS_SUPERDROPS_THERMODYNAMIC_EQUATIONS_HPP_
