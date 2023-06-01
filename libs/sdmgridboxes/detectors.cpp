// Author: Clara Bayley
// File: "detectors.cpp"
/* Fucntionality of detectors
(e.g. of SDM processes) in gridboxes
which copy data from detections
into 'logbooks' */

#include "./detectors.hpp"

double AccumPrecipDetector::
    accumulated_precipitation(const Superdrop &drop) const
{
  std::cout << "PRECIP!\n";
  return 0.0;
}

std::shared_ptr<Detectors> DetectorsInstallation::
    install_precipitation_detectors(const std::shared_ptr<Detectors> detectors,
                                    const unsigned int gbxindex,
                                    const Maps4GridBoxes &gbxmaps) const
/* if upper z boundary of gbx is <= precip_zlim install
a detector to detect accumulated precipitation */
{
  constexpr double precip_zlim(50 / dlc::COORD0); // (dimless) maximum z coord of gbxs that detect precipitation

  if (gbxmaps.get_bounds_z(gbxindex).second <= precip_zlim)
  {
    detectors -> install_accumprecip_detector(logbooks, gbxindex);
  }

  return detectors;
}

std::shared_ptr<Detectors> DetectorsInstallation::
operator()(const unsigned int gbxindex,
           const Maps4GridBoxes &gbxmaps) const
/* operator creates a unique pointer to a
detectors struct and installs certain
types of detector in it */
{
  auto detectors = std::make_shared<Detectors>();

  detectors = install_precipitation_detectors(detectors,
                                              gbxindex,
                                              gbxmaps);

  return detectors;
}