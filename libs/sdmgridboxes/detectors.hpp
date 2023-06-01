// Author: Clara Bayley
// File: "detectors.hpp"
/* Header file for functions and
structures related to detectors
(e.g. of SDM processes) in gridboxes
which copy data from detections
into 'logbooks' */

#ifndef DETECTORS_HPP
#define DETECTORS_HPP

#include <memory>

#include "./logbooks.hpp"
#include "./maps4gridboxes.hpp"
#include "../claras_SDconstants.hpp"
#include "../superdrop_solver/superdrop.hpp"

namespace dlc = dimless_constants;

struct AccumPrecipDetector
/* detector which stores the value of
accumulated precipitation in an entry of
a logbook controlled by the
EntryInLogbook instance */
{
private:
  EntryInLogbook<double> manage_entry;

  double precipitation(const Superdrop &drop) const;

public:
  KOKKOS_INLINE_FUNCTION ~AccumPrecipDetector() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  AccumPrecipDetector() : manage_entry() {} // also a Kokkos requirement for a (dual)View
  /* initialise without a logbook */

  KOKKOS_INLINE_FUNCTION
  AccumPrecipDetector(const std::shared_ptr<Logbook<double>> logbook,
                      const unsigned int gbxindex)
      : manage_entry(logbook, gbxindex) {}
  /* initialise manage_entry with a
  logbook with tag 'gbxindex'*/

  void operator()(const Superdrop &drop) const
  /* if detector has a logbook, use manage_entry to
  store accumlated precipitation in it */
  {
    if (manage_entry.get_logbook())
    {
      manage_entry.increment_by(precipitation(drop));
    }
  }
};

class Detectors
/* Detectors stores various detector types and 
a reference to logbook instances found in
'DetectionLogbooks'. Detectors calls default 
constructor of detector types upon construction.
Modification of a detector can then
be done by calling the appropriate install_[...] 
function. Likewise detcetor can be used through 
appropriate detect_[...] function */
{
private:
  AccumPrecipDetector accpp_dtr;

public:
  KOKKOS_INLINE_FUNCTION Detectors() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Detectors() = default; // Kokkos requirement for a (dual)View

  void install_accumprecip_detector(
      const std::shared_ptr<Logbook<double>> accpp_logbook,
      const unsigned int gbxindex)
  /* install accumulated precipitation detector
  (by instanting detector with an entry in the
  accpp logbook that has tag 'gbxindex') */
  {
    accpp_dtr = AccumPrecipDetector(accpp_logbook, gbxindex);
  }

  void detect_precipitation(const Superdrop &drop) const
  /* use operators of precipitation detectors
  to detect precipitation  */
  {
    accpp_dtr(drop); 
  }
};

struct DetectorsInstallation
/* DetectorsInstallation operator()
returns a smart pointer to a detectors
instance that may modify data in
vectors pointed to by logbooks */
{
private:
  const DetectorLogbooks &logbooks;
  const Maps4GridBoxes &gbxmaps;

  std::shared_ptr<Detectors>
  install_precipitation_detectors(const std::shared_ptr<Detectors> detectors,
                                  const unsigned int gbxindex) const;
  /* if upper z boundary of gbx is <= precip_zlim install
  a detector to detect accumulated precipitation */

  std::shared_ptr<Detectors>
  install_detectors(std::shared_ptr<Detectors> detectors,
                    const unsigned int gbxindex) const;
  /* operator installs certain types of detector in
  detectors struct given its pointer */

public:
  DetectorsInstallation(const DetectorLogbooks &logbooks,
                        const Maps4GridBoxes &gbxmaps)
      : logbooks(logbooks), gbxmaps(gbxmaps) {}

  std::shared_ptr<Detectors> operator()(const unsigned int gbxindex) const
  /* operator creates a unique pointer to a detectors struct
  and installs certain types of detector in it according
  to install_detectors function */
  {
    auto detectors = std::make_shared<Detectors>();
    
    return install_detectors(detectors, gbxindex);
  }
};

#endif // DETECTORS_HPP