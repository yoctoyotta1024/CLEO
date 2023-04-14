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

#include "superdrop_solver/superdrop.hpp"
#include "./singlevarstorage.hpp"

struct SDMomentsStorage
{
  const double scalefac; // scale factor to convert dimensionless masses to grams
 
  TwoDStorage<double> massmoment0zarr;
  TwoDStorage<double> massmoment1zarr;
  TwoDStorage<double> massmoment2zarr;

  SDMomentsStorage(FSStore &store, const unsigned int maxcsize,
                   const unsigned int ngridboxes)
      : scalefac(pow(dlc::R0, 3.0) * dlc::RHO0 * 1000),
        massmoment0zarr(store, maxcsize, "massmoment0", "<f8",
                        " ", 1, ngridboxes),
        massmoment1zarr(store, maxcsize, "massmoment1", "<f8",
                        "g", scalefac, ngridboxes),
        massmoment2zarr(store, maxcsize, "massmoment2", "<f8",
                        "g^2", pow(scalefac, 2.0), ngridboxes){};
};

double massnthmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx,
                    const double nth_moment);
/* calculates the nth moment of the (real) droplet mass distirbution
given by the superdrops in the span passed as an argument */

#endif // SDMOMENTSSTORAGE_HPP 