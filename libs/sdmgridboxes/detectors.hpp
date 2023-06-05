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

struct SurfPrecipDetector
/* detector which stores the value of
accumulated precipitation over a time 
duration in an entry of a logbook controlled
by the EntryInLogbook instance */
{
private:
  EntryInLogbook<double> manage_entry;

  double precip_mass(const Superdrop &drop) const
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

public:
  KOKKOS_INLINE_FUNCTION ~SurfPrecipDetector() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  SurfPrecipDetector() : manage_entry() {} // also a Kokkos requirement for a (dual)View
  /* initialise without a logbook */

  KOKKOS_INLINE_FUNCTION
  SurfPrecipDetector(const std::shared_ptr<Logbook<double>> logbook,
                      const unsigned int gbxindex)
      : manage_entry(logbook, gbxindex) {}
  /* initialise manage_entry with a
  logbook with tag 'gbxindex'*/

  void operator()(const Superdrop &drop) const
  /* if detector has a logbook, use manage_entry to store
  accumlated precipitation over some duration in it */
  {
    if (manage_entry.get_logbook())
    {
      manage_entry.increment_by(precip_mass(drop));
    }
  }
};

class Detectors
/* Detectors stores various detector types and 
calls default constructor of them upon construction.
Installation/Modification of a detector can then
be done by calling the appropriate install_[...] 
function. Likewise detcetor can be used through 
appropriate detect_[...] function */
{
private:
  SurfPrecipDetector detect_surfprecip;

public:
  KOKKOS_INLINE_FUNCTION Detectors() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Detectors() = default; // Kokkos requirement for a (dual)View

  void install_surfprecip_detector(
      const std::shared_ptr<Logbook<double>> surfpp_logbook,
      const unsigned int gbxindex)
  /* install accumulated precipitation detector
  (by instanting detector with an entry in the
  surfprecip logbook that has tag 'gbxindex') */
  {
    detect_surfprecip = SurfPrecipDetector(surfpp_logbook, gbxindex);
  }

  void detect_precipitation(const Superdrop &drop) const
  /* use operators of precipitation detectors
  to detect precipitation  */
  {
    detect_surfprecip(drop); 
  }
};

#endif // DETECTORS_HPP