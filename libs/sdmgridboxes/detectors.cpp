// Author: Clara Bayley
// File: "detectors.cpp"
/* Fucntionality of detectors
(e.g. of SDM processes) in gridboxes */

#include "./detectors.hpp"

void Detectors::precipitation(const unsigned int gbxindex, const Superdrop drop)
{
  if (unique_pointer_to_position_in_accumulatedprecip) // <- if detector points to a location (e.g. in accumulatde precip vector)
  {
    if (drop.coord3 < 0.0)
    {
      position_in_accumulatedprecip += drop.mass();
      // accum_precip_radii.pushback(drop.radius()); // for observing raindrop distribution
      // accum_precip_ndrops += 1;
    }

  }
  
}