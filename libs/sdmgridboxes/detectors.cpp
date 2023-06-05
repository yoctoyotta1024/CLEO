// Author: Clara Bayley
// File: "detectors.cpp"
/* functions related to detectors
(e.g. of SDM processes) */

#include "./detectors.hpp"

double SurfPrecipDetector::precip_mm(const double gbxarea,
                                     const Superdrop &drop) const
/* returns (dimless) amount of precipitation calculated from
volume of liquid from (real) droplets when superdroplet
coord3 is below 0.0. Amount of precip in millimeters
= return_value * R0^3 / COORD0^2 * 1000 */
{
  if (drop.coord3 < 0.0)
  {
    return drop.vol_liq() * drop.eps / gbxarea; // amount of precip (dimless)
  }

  return 0.0;
}