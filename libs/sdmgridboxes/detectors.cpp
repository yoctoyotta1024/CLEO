// Author: Clara Bayley
// File: "detectors.cpp"
/* Fucntionality of detectors
(e.g. of SDM processes) in gridboxes
which copy data from detections
into 'logbooks' */

#include "./detectors.hpp"

double AccumPrecipDetector::
    precip_mass(const Superdrop &drop) const
/* returns (dimless) mass of precipitation
calulated as mass of (real) droplets
when superdroplet is below coord3 = 0.0 */
{
  if (drop.coord3 < 0.0)
  {
    return drop.mass() * drop.eps;
  }
  
  return 0.0;
}