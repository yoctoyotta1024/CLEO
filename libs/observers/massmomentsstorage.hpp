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

struct MassMom012Storages
/* 2D zarr stores for 0th, 1st and 2nd mass moments */
{
private:
  const double sf = pow(dlc::R0, 3.0) * dlc::RHO0 * 1000; // scale factor to convert dimensionless masses to grams

public:
  TwoDStorage<double> mom0zarr;
  TwoDStorage<double> mom1zarr;
  TwoDStorage<double> mom2zarr;

  MassMom012Storages(FSStore &store, const unsigned int maxchunk,
                       const unsigned int ngridboxes,
                       const std::string mom0name,
                       const std::string mom1name,
                       const std::string mom2name)
      : mom0zarr(store, maxchunk, mom0name,
                 "<f8", " ", 1.0, ngridboxes),
        mom1zarr(store, maxchunk, mom1name,
                 "<f8", "g", sf, ngridboxes),
        mom2zarr(store, maxchunk, mom2name,
                 "<f8", "g^2", pow(sf, 2.0), ngridboxes){};
};

struct MassMomStorages : MassMom012Storages
{
  MassMomStorages(FSStore &store, const unsigned int maxchunk,
           const unsigned int ngridboxes)
      : MassMom012Storages(store, maxchunk, ngridboxes,
                             "mom0", "mom1", "mom2"){};
};

struct RainMassMomStorages : MassMom012Storages
{
  RainMassMomStorages(FSStore &store, const unsigned int maxchunk,
               const unsigned int ngridboxes)
      : MassMom012Storages(store, maxchunk, ngridboxes,
                             "rainmom0", "rainmom1", "rainmom2"){};
};

double massmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx,
                  const double nth_moment);
/* calculates the nth moment of the (real) droplet mass distirbution
given by all the superdrops in the span passed as an argument */

double rainmassmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx,
                      const double nth_moment);
/* calculates the nth moment of the (real) raindroplet
mass distirbution given by all the superdrops which
have radius >= rlim in the span passed as an argument */

double surface_precipitation(const GridBox &gbx, const double coord3lim);
/* calculates mm of precipitation in a gridbox
from mass of all superdrops which have
radius >= rlim and coord3 <= zlim  */

#endif // MASSMOMENTSSTORAGE_HPP