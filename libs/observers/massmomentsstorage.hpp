// Author: Clara Bayley
// File: "massmomentsstorage.hpp"
/* structs/classes to create an observer for the moments of
the superdroplet mass distribution that writes
into 1 dimensional array(s)
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#ifndef MASSMOMENTSSTORAGE_HPP
#define MASSMOMENTSSTORAGE_HPP

#include <span>
#include <cmath>

#include "../claras_SDconstants.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "./singlevarstorage.hpp"

namespace dlc = dimless_constants;

struct MassMomentsStorage
{
  const double sf; // scale factor to convert dimensionless masses to grams

  TwoDStorage<double> massmom0zarr;
  TwoDStorage<double> massmom1zarr;
  TwoDStorage<double> massmom2zarr;

  MassMomentsStorage(FSStore &store, const unsigned int maxchunk,
                   const unsigned int ngridboxes)
      : sf(pow(dlc::R0, 3.0) * dlc::RHO0 * 1000),
        massmom0zarr(store, maxchunk,
                     "massmom0", "<f8", " ", 1.0, ngridboxes),
        massmom1zarr(store, maxchunk,
                     "massmom1", "<f8", "g", sf, ngridboxes),
        massmom2zarr(store, maxchunk,
                     "massmom2", "<f8", "g^2", pow(sf, 2.0), ngridboxes){};
};

double massmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx,
                     const double nth_moment);
/* calculates the nth moment of the (real) droplet mass distirbution
given by all the superdrops in the span passed as an argument */

#endif // MASSMOMENTSSTORAGE_HPP