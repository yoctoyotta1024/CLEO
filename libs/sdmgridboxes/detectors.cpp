// Author: Clara Bayley
// File: "detectors.cpp"
/* Fucntionality of detectors
(e.g. of SDM processes) in gridboxes
which copy data from detections
into 'logbooks' */

#include "./detectors.hpp"

double AccumPrecipDetector::
    precipitation(const Superdrop &drop) const
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

std::shared_ptr<Detectors> DetectorsInstallation::
    install_precipitation_detectors(const std::shared_ptr<Detectors> detectors,
                                    const unsigned int gbxindex) const
/* if upper z boundary of gbx is <= precip_zlim install
a detector to detect accumulated precipitation */
{
  constexpr double precip_zlim(50 / dlc::COORD0); // (dimless) maximum z coord of gbxs that detect precipitation

  if (gbxmaps.get_bounds_z(gbxindex).second <= precip_zlim)
  {
    detectors->install_accumprecip_detector(logbooks.accpp, gbxindex);
  }
  return detectors;
}

std::shared_ptr<Detectors> DetectorsInstallation::
    install_detectors(std::shared_ptr<Detectors> detectors,
                      const unsigned int gbxindex) const
/* operator installs certain types of detector in
detectors struct given its pointer */
{
  detectors = install_precipitation_detectors(detectors, gbxindex);

  return detectors;
}