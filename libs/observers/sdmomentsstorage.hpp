// Author: Clara Bayley
// File: "sdmomentsstorage.hpp"
/* structs/classes to create an observer for the moments of
the superdroplet mass distribution that writes
into 1 dimensional array(s)
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#ifndef SDMOMENTSSTORAGE_HPP
#define SDMOMENTSSTORAGE_HPP

#include <span>
#include <cmath>

#include "../claras_SDconstants.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "./singlevarstorage.hpp"

namespace dlc = dimless_constants;

struct SDMomentsStorage
{
  const double sf; // scale factor to convert dimensionless masses to grams

  TwoDStorage<double> massmom0zarr;
  TwoDStorage<double> massmom1zarr;
  TwoDStorage<double> massmom2zarr;

  SDMomentsStorage(FSStore &store, const unsigned int maxchunk,
                   const unsigned int ngridboxes)
      : sf(pow(dlc::R0, 3.0) * dlc::RHO0 * 1000),
        massmom0zarr(store, maxchunk,
                     "massmom0", "<f8", " ", 1.0, ngridboxes),
        massmom1zarr(store, maxchunk,
                     "massmom1", "<f8", "g", sf, ngridboxes),
        massmom2zarr(store, maxchunk,
                     "massmom2", "<f8", "g^2", pow(sf, 2.0), ngridboxes){};
};

double massnthmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx,
                     const double nth_moment);
/* calculates the nth moment of the (real) droplet mass distirbution
given by the superdrops in the span passed as an argument */

#endif // SDMOMENTSSTORAGE_HPP