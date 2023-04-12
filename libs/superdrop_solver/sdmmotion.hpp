// Author: Clara Bayley
// File: "sdmmotion.hpp"
/* Header file for functions related to
moving superdroplets (updating their
coordinates) */

#ifndef SDMMOTION_HPP
#define SDMMOTION_HPP

#include <span>

#include "./superdrop.hpp"

struct SdmMotion
{
  SdmMotion(){};

  void move_superdroplets(std::span<SuperdropWithGbxindex> span4SDsinGBx,
                 const double w, const double u, const double v) const;
};

#endif // SDMMOTION_HPP