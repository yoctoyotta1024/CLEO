// Author: Clara Bayley
// File: "detectors_ptr.cpp"
/* Functionality for some types
obeying the CreateDetectorsPtr concept */

#include "./detectors_ptr.hpp"

std::shared_ptr<Detectors> PrecipDetectorsPtr::
    install_precip_detectors(const std::shared_ptr<Detectors> detectors,
                                    const DetectorLogbooks &logbooks,
                                    const unsigned int gbxindex) const
/* if upper z boundary of gbx is <= precip_zlim install a detector
to detect accumulated precipitation at surface over 1 timestep */
{
  // constexpr double precip_zlim(50 / dlc::COORD0); // (dimless) maximum z coord of gbxs that detect precipitation
  const double precip_zlim(gbxmaps.get_bounds_z(0).second); // z boundary of lowest layer gbxs

  if (gbxmaps.get_bounds_z(gbxindex).second <= precip_zlim)
  {
    detectors->install_surfprecip_detector(logbooks.surfpp, gbxindex);
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