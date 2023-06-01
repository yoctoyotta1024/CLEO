// Author: Clara Bayley
// File: "detectors.hpp"
/* Header file for functions and
structures related to detectors
(e.g. of SDM processes) in gridboxes */

#ifndef DETECTORS_HPP
#define DETECTORS_HPP

#include <vector>

#include "superdrop_solver/superdrop.hpp"


struct DetectionRecords
{
  std::vector<double> accumulatedprecip; // <- each value here is accum precip of a gbx
};

struct Detectors
{
  // unique_pointer_to_position_in_accumulatedprecip;

  void precipitation(const unsigned int gbxindex, const Superdrop drop);
};


#endif // DETECTORS_HPP