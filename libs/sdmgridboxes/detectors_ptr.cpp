// Author: Clara Bayley
// File: "detectors_ptr.cpp"
/* Functionality for some types
obeying the CreateDetectorsPtr concept */

#include "./detectors_ptr.hpp"

std::shared_ptr<Detectors> PrecipDetectorsPtr::
    install_precip_detectors(const std::shared_ptr<Detectors> detectors,
                                    const DetectorLogbooks &logbooks,
                                    const unsigned int gbxindex) const
/* if upper z boundary of gbx is <= precip_zlim install
a detector to detect accumulated precipitation */
{
  constexpr double precip_zlim(500 / dlc::COORD0); // (dimless) maximum z coord of gbxs that detect precipitation

  if (gbxmaps.get_bounds_z(gbxindex).second <= precip_zlim)
  {
    detectors->install_accumprecip_detector(logbooks.accumprecip, gbxindex);
  }
  return detectors;
}

std::shared_ptr<Detectors> PrecipDetectorsPtr::
    install_detectors(std::shared_ptr<Detectors> detectors,
                      const DetectorLogbooks &logbooks,
                      const unsigned int gbxindex) const
/* operator installs certain types of detector in
detectors struct given its pointer */
{
  detectors = install_precip_detectors(detectors, logbooks, gbxindex);

  return detectors;
}