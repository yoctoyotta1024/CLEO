// Author: Clara Bayley
// File: "sdmmotion.hpp"
/* Header file for functions related to
moving superdroplets (updating their
coordinates) */

#ifndef SDMMOTION_HPP
#define SDMMOTION_HPP

#include "./superdrop.hpp"

struct SdmMotion
{
  SdmMotion(){};

  void move_superdroplet(const ThermoState &state,
                         Superdrop &superdrop) const;
};

#endif // SDMMOTION_HPP