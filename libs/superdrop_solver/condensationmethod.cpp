// Author: Clara Bayley
// File: condensationmethod.cpp
/* Functionality for modelling condensation-
diffusional growth/shrinking of superdroplets. Equations
referenced as (eqn [X.YY]) are from "An Introduction
To Clouds From The Microscale to Climate" by
Lohmann, Luond and Mahrt, 1st edition. */

#include "condensationmethod.hpp"

std::pair<double, double> CondensationMethod::diffusion_factors(const double press,
                                                                const double temp,
                                                                const double psat) const
/* Calculate dimensionless Fkl and Fdl
  heat and vapour diffusion factors in
  equation for radial growth of droplets
  according to equations from "An Introduction
  To Clouds...." (see note at top of file).
  fkl is first item of returned pair, fdl is second. */
{
  constexpr double A = 7.11756e-5;                            // coefficient for T^2 in T*[eq.7.24]
  constexpr double B = 4.38127686e-3;                         // coefficient for T in T*[eq.7.24]
  constexpr double LATENT_RGAS_V = DC::LATENT_V / DC::RGAS_V; // for fkl diffusion factor calc
  constexpr double D = 4.012182971e-5;                        // constants in equation [eq.7.26]

  const double TEMP = temp * dlc::TEMP0;
  const double PRESS = press * dlc::P0;
  const double PSAT = psat * dlc::P0;

  const double THERMK = A * pow(TEMP, 2.0) + TEMP * B;                 // K*TEMP with K from [eq.7.24] (for fkl)
  const double DIFFUSE_V = (D / PRESS * pow(TEMP, 1.94)) / DC::RGAS_V; // 1/R_v * D_v from [eq 7.26] (for fdl)

  const double fkl = (LATENT_RGAS_V / TEMP - 1.0) * DC::LATENT_V / (THERMK * dlc::F0); // fkl eqn [7.23]
  const double fdl = TEMP / (DIFFUSE_V * PSAT) / dlc::F0;                              // fdl eqn [7.25]

  return std::pair<double, double>(fkl, fdl);
}

double CondensationMethod::superdroplet_growth_by_condensation(const double press,
                                                               const double temp, const double psat,
                                                               const double s_ratio, const double delt,
                                                               const ImplicitEuler &impliciteuler, Superdrop &drop) const
/* update superdroplet radius due to radial growth/shrink
  via condensation and diffusion of water vapour according
  to equations from "An Introduction To Clouds...." (see
  note at top of file). Then return mass of liquid that
  condensed/evaporated onto droplet. New radius is calculated
  using impliciteuler method which iterates condensation-diffusion
  ODE given the previous radius. */
{
  /* n.b. Structured bindings require C++17. Can use std::tie for C++11 */
  const double dmdt_const = 4.0 * M_PI * drop.get_solute()->rho_l * pow(dlc::R0, 3.0);
  const double akoh = drop.akohler_factor(temp);
  const double bkoh = drop.bkohler_factor();
  const auto [fkl, fdl] = diffusion_factors(press, temp, psat);

  /* do not pass r by reference here!! copy value into iterator */
  double newradius = impliciteuler.implicitmethod_forcondensation(s_ratio,
                                                                  akoh, bkoh, fkl,
                                                                  fdl, drop.radius); // timestepping eqn [7.28] forward
  double delta_radius = drop.change_radius(newradius);
  double mass_condensed = (dmdt_const * pow(drop.radius, 2.0) * drop.eps * delta_radius); // eqn [7.22] * delta t

  return mass_condensed;
}

void CondensationMethod::condensation_alters_thermostate(ThermoState &state,
                                                          const double tot_rho_condensed) const
/* change the thermodynamic variables (temp, qv and qc) of
ThermoState state given the total change in condensed
water mass per volume during time interval delt */
{
  const double delta_qcond = tot_rho_condensed / dlc::Rho_dry;
  const double delta_qvap = -(delta_qcond);
  const double delta_temp = (dlc::Latent_v /
                             moist_specifc_heat(state.qvap, state.qcond)) *
                            delta_qcond;

  state.temp += delta_temp;
  state.qvap += delta_qvap;
  state.qcond += delta_qcond;
}

void CondensationMethod::condensation_onto_superdroplets(std::span<SuperdropWithGbxindex> span4SDsinGBx,
                                                         ThermoState &state) const
/* Change to superdroplet radii and temp, qv and
qc due to sum of radii changes via diffusion and
condensation of water vapour during timestep delt.
Using equations from "An Introduction To
Clouds...." (see note at top of file) */
{
  const double psat = saturation_pressure(state.temp);
  const double s_ratio = supersaturation_ratio(state.press, state.qvap, psat);

  double tot_rho_condensed = 0.0; // cumulative change to liquid mass in parcel volume

  /* superdroplet radii changes */
  for (auto &SDinGBx : span4SDsinGBx)
  {
    const double delta_mass_condensed = superdroplet_growth_by_condensation(state.press, state.temp,
                                                                            psat, s_ratio, delt,
                                                                            impliciteuler, SDinGBx.superdrop);
    const double VOLUME = state.volume * pow(dlc::COORD0, 3.0); // volume in which condensation occurs [m^3]
    tot_rho_condensed += (delta_mass_condensed / VOLUME);       // drho_condensed_vapour/dt * delta t
  }

  /* resultant effect on thermodynamic state */
  if (doAlterThermo)
  {
    condensation_alters_thermostate(state, tot_rho_condensed);
  }
}