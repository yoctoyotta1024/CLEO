// Author: Clara Bayley
// File: condensationmethod.cpp
/* Class containing method for condensation-
diffusional growth/shrinking of superdroplets.
Equations referenced as (eqn [X.YY]) are
from "An Introduction To Clouds From The
Microscale to Climate" by Lohmann, Luond
and Mahrt, 1st edition. */

#ifndef CONDENSATIONMETHOD_HPP
#define CONDENSATIONMETHOD_HPP

#include <random>
#include <span>
#include <functional>
#include <concepts>

#include "../claras_SDconstants.hpp"
#include "./thermodynamic_equations.hpp"
#include "./impliciteuler.hpp"
#include "./thermostate.hpp"
#include "./superdrop.hpp"
#include "./sdmprocess.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

class CondensationMethod
/* class for method to grow/shrink superdroplets due to
condensation/evaporation during some constant time interval.
Optionally also implements the resultant thermodynamic
changes to a (ThermoState) thermodynamic state */
{
private:
  const bool doAlterThermo;               // whether to make condensation alter ThermoState or not
  const double delt;                 // dimensionless time interval during which condenstaion occurs
  const ImplicitEuler impliciteuler; // method to integrate condensation equation

  std::pair<double, double> diffusion_factors(const double press, const double temp,
                                              const double psat) const;
  /* Calculate dimensionless Fkl and Fdl
  heat and vapour diffusion factors in
  equation for radial growth of droplets
  according to equations from "An Introduction
  To Clouds...." (see note at top of file).
  fkl is first item of returned pair, fdl is second. */

  double superdroplet_growth_by_condensation(const double press,
                                             const double temp, const double psat,
                                             const double s_ratio, const double delt,
                                             const ImplicitEuler &impliciteuler, Superdrop &drop) const;
  /* update superdroplet radius due to radial growth/shrink
  via condensation and diffusion of water vapour according
  to equations from "An Introduction To Clouds...." (see
  note at top of file). Then return mass of liquid that
  condensed/evaporated onto droplet. New radius is calculated
  using impliciteuler method which iterates condensation-diffusion
  ODE given the previous radius. */

  void condensation_alters_thermostate(ThermoState &state,
                                        const double tot_rho_condensed) const;
  /* change the thermodynamic variables (temp, qv and qc)
  given the total change in condensed mass per
  parcel volume during timestep delt */

  void condensation_onto_superdroplets(std::span<SuperdropWithGbxindex> span4SDsinGBx,
                                       ThermoState &state) const;
  /* Change to superdroplet radii and temp, qv and
  qc due to sum of radii changes via diffusion and
  condensation of water vapour during timestep delt.
  Using equations from "An Introduction To
  Clouds...." (see note at top of file) */

public:
  CondensationMethod(const bool doAlterThermo, const double delt,
                     const ImplicitEuler impliciteuler)
      : doAlterThermo(doAlterThermo),
        delt(delt),
        impliciteuler(impliciteuler) {}

  CondensationMethod(const bool doAlterThermo, const double delt,
                     const unsigned maxiters, const double rtol,
                     const double atol)
      : doAlterThermo(doAlterThermo),
        delt(delt),
        impliciteuler(maxiters, delt, rtol, atol) {}

  inline void operator()(const int currenttimestep,
                         std::span<SuperdropWithGbxindex> span4SDsinGBx,
                         ThermoState &state,
                         std::mt19937 &gen) const
  /* this operator is used as an "adaptor" for using a run_step
  function in order to call condensation_onto_superdroplets.
  (*hint* run_step is usually found within a type that
  satisfies the SdmProcess concept) */
  {
    condensation_onto_superdroplets(span4SDsinGBx, state);
  }
};

SdmProcess auto CondensationProcess(const int interval,
                                    const std::function<double(int)> int2time,
                                    const bool doAlterThermo,
                                    const unsigned maxiters,
                                    const double rtol,
                                    const double atol)
/* constructs SdmProcess for condensation with constant timestep 'interval'
given a function to convert the interval to a (dimensionless) time
and the arguments required to construct the condensation method */
{
  const double dimlesststep = int2time(interval);
  return ConstTstepProcess{interval, CondensationMethod(doAlterThermo, dimlesststep,
                                                        maxiters, rtol, atol)};
}

#endif // CONDENSATIONMETHOD_HPP