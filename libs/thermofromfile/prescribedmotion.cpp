// Author: Clara Bayley
// File: "prescribedmotion.cpp"
/* Functionality to update
superdroplet positions 
according to a prescribed 2D flow */

#include "prescribedmotion.hpp"

Prescribed2DFlow::Prescribed2DFlow(const double zlength,
                                   const double xlength,
                                   const double wmax,
                                   const std::function<double(ThermoState)> rhotilda)
    : ztilda(zlength / std::numbers::pi),         // 1/wavenumber given dimensionless wavelength
      xtilda(xlength / (2.0 * std::numbers::pi)), // 1/wavenumber given dimensionless wavelength
      wamp(2.0 * wmax),                           // amplitude of velocioty variations
      rhotilda(rhotilda)                          // normaliseddry air density
{
}

double Prescribed2DFlow::prescribed_wvel(const ThermoState &state,
                                         const double zcoord, const double xcoord) const
{
  return wamp / rhotilda(state) *
         std::sin(zcoord / ztilda) * std::sin(xcoord / xtilda);
}

double Prescribed2DFlow::prescribed_uvel(const ThermoState &state,
                                         const double zcoord, const double xcoord) const
{
  return wamp / rhotilda(state) * xtilda / ztilda *
         std::cos(zcoord / ztilda) * std::cos(xcoord / xtilda);
}

void MoveWith2DPrescribedFlow::
    change_superdroplet_coords(const Maps4GridBoxes &gbxmaps,
                               const GridBox &gbx,
                               Superdrop &drop) const
/* Use predictor-corrector scheme from Grabowksi et al. 2018
(similar to Arabas et al. 2015) to update a SD position.
The velocity required for this scheme is determined
from the PrescribedFlow2D instance */
{
  auto deltas = predictor_corrector(gbx.state, drop.coord3, drop.coord1);
  // auto deltas = leapfrog(gbx.state, drop.coord3, drop.coord1);
  
  cfl_criteria(gbxmaps, gbx.gbxindex, deltas.first, deltas.second, 0.0);

  drop.coord3 += deltas.first;
  drop.coord1 += deltas.second; 
}

std::pair<double, double>
MoveWith2DPrescribedFlow::predictor_corrector(const ThermoState &state,
                                         const double coord3,
                                         const double coord1) const
/* returns change in (z,x) coordinates = (delta3, delta1)
obtained using predictor-corrector method and velocities
calculated from a Prescribed2DFlow */
{ 
  const double vel3 = flow2d.prescribed_wvel(state, coord3, coord1); // w wind from prescribed 2D flow
  const double vel1 = flow2d.prescribed_uvel(state, coord3, coord1); // u wind from prescribed 2D flow

  const double pred3(coord3 + vel3 * delt);
  const double pred1(coord1 + vel1 * delt);

  const double corrvel3 = flow2d.prescribed_wvel(state, pred3, pred1); // w wind from prescribed 2D flow
  const double corrvel1 = flow2d.prescribed_uvel(state, pred3, pred1); // w wind from prescribed 2D flow

  const double delta3((vel3 + corrvel3) * (delt / 2));
  const double delta1((vel1 + corrvel1) * (delt / 2));

  return std::pair<double, double>(delta3, delta1);
}

std::pair<double, double>
MoveWith2DPrescribedFlow::leapfrog(const ThermoState &state,
                              const double coord3,
                              const double coord1) const
/* returns change in (z,x) coordinates = (delta3, delta1)
obtained using a simple leapfrog method and velocities
calculated from a Prescribed2DFlow */
{ 
  const double vel1 = flow2d.prescribed_uvel(state, coord3, coord1); // w wind from prescribed 2D flow
  const double pred1(coord1 + vel1 * delt);

  const double vel3 = flow2d.prescribed_wvel(state, coord3, pred1); // u wind from prescribed 2D flow
  
  const double delta3(vel3 * delt);
  const double delta1(vel1 * delt);

  return std::pair<double, double>(delta3, delta1);
}