/*
 * ----- CLEO -----
 * File: cleoconstants.hpp
 * Project: libs
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 14th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * header file with namespaces for (physical) constants used by CLEO
 * Note: All letters in CAPITALS indicates constants have dimensions
 * (or are naturally dimensionless as for IONIC). Only first letter
 * capitalised indicates value is dimensionless (version of a) constant
 */

#ifndef CLEOCONSTANTS_HPP 
#define CLEOCONSTANTS_HPP 

#include <limits>

namespace dimmed_constants
/* namespace containing values of
constants with dimensions */
{
  constexpr double G = 9.80665;                   // acceleration due to gravity [m/s^2]
  constexpr double RGAS_UNIV = 8.314462618;       // universal molar gas constant [J/Kg/K]
  constexpr double MR_WATER = 0.01801528;         // molecular mass of water [Kg/mol]
  constexpr double MR_DRY = 0.028966216;          // molecular mass of dry air [Kg/mol]
  constexpr double RGAS_DRY = RGAS_UNIV / MR_DRY; // specific gas constant for dry air [J/Kg/K]
  constexpr double RGAS_V = RGAS_UNIV / MR_WATER; // specific gas constant for water [J/Kg/K]

  constexpr double LATENT_V = 2500930; // specific latent heat of vapourisation of water [J/Kg]  (IAPWS97 at 273.15K)
  constexpr double CP_DRY = 1004.64;   // specific heat capacity (dry) air at constant pressure [J/Kg/K] ~1.400*cv_dry (ICON)
  constexpr double CP_V = 1865.01;     // specific heat capacity of water vapour [J/Kg/K] (IAPWS97 at 273.15K)
  constexpr double C_L = 4192.664;     // specific heat capacity of liquid water[J/Kg/K] (ICON c_l = (3.1733 + 1.0) * cp_dry)

  constexpr double RHO_DRY = 1.177;        // density of dry air [Kg/m^3] (at 300K)
  constexpr double RHO_L = 998.203;        // density of liquid water condensing [kg/m^3] (water at 293K from SCALE-SDM)
  constexpr double DYNVISC = 18.45 * 1e-6; // dynamic viscosity of air [Pa s] (at 300K)

  constexpr double RHO_SOL = 2016.5;    // density of (dry) areosol [Kg/m^3] (NaCl from SCALE-SDM)
  constexpr double MR_SOL = 0.05844277; //  molecular mass of areosol [Kg/mol] (NaCl=0.058 from SCALE-SDM)
  constexpr int IONIC = 2;            //  degree ionic dissociation (van't Hoff factor) [dimensionless]
}

namespace dimless_constants
{
  /* constants for using characterstic scales of time, velocity,
  temperature etc. (TIME0, TEMP0, P0 etc.) in order to make variables
  dimensionless. Namespace also includes the dimensionless equivalents
  of some members of the dimmed_constants namespace*/
  namespace DC = dimmed_constants;

  /* characterstic scales */
  constexpr double W0 = 1.0;            // characteristic velocity [m/s]
  constexpr double TIME0 = 1000.0;      // timescale [s]
  constexpr double COORD0 = TIME0 * W0; // coordinate grid scale [m]

  constexpr double CP0 = DC::CP_DRY;              // Heat capacity [J/Kg/K]
  constexpr double MR0 = DC::MR_DRY;              // molecular molar mass [Kg/mol]
  constexpr double R0 = 1e-6;                     // droplet radius lengthscale [m]

  constexpr double P0 = 100000.0;                    // pressure [Pa]
  constexpr double TEMP0 = 273.15;                   // temperature [K]
  constexpr double RHO0 = P0 / (CP0 * TEMP0);        // density [Kg/m^3]
  constexpr double MASS0 = R0 * R0 * R0 * RHO0;      // mass [Kg]
  constexpr double MASS0grams = MASS0 * 1000;        // mass [g]
  constexpr double F0 = TIME0 / (RHO0 * R0 * R0);    // droplet condensation-diffusion factors []

  /* dimensionaless constants */
  constexpr double Mr_ratio = DC::MR_WATER / DC::MR_DRY;
  constexpr double Cp_dry = DC::CP_DRY / CP0;
  constexpr double Cp_v = DC::CP_V / CP0;
  constexpr double C_l = DC::C_L / CP0;
  constexpr double Latent_v = DC::LATENT_V / (TEMP0 * CP0);
  constexpr double Rgas_dry = DC::RGAS_DRY / CP0;
  constexpr double Rgas_v = DC::RGAS_V / CP0;
  constexpr double Rho_dry = DC::RHO_DRY / RHO0;
  constexpr double Rho_l = DC::RHO_L / RHO0;
  constexpr double Rho_sol = DC::RHO_SOL / RHO0;
  constexpr double Mr_sol = DC::MR_SOL / MR0;
  constexpr int IONIC = DC::IONIC;
  
  /* Domain, Gridbox and/or Superdroplet constants */
  constexpr unsigned int OUTOFDOMAIN = std::numeric_limits<unsigned int>::max(); // value for gbxindex that indicates superdrop is out of domain boundaries (ie. no longer exists)
}

#endif // CLEOCONSTANTS_HPP