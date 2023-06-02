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

#include "./singlevarstorage.hpp"
#include "../claras_SDconstants.hpp"
#include "superdrop_solver/superdrop.hpp"

namespace dlc = dimless_constants;

struct MassMom012Storages
/* 2D zarr stores for 0th, 1st and 2nd mass moments */
{
  TwoDStorage<double> mom0zarr;
  TwoDStorage<double> mom1zarr;
  TwoDStorage<double> mom2zarr;

  MassMom012Storages(FSStore &store, const unsigned int maxchunk,
                     const unsigned int ngbxs,
                     const std::string name0,
                     const std::string name1,
                     const std::string name2)
      : mom0zarr(store, maxchunk, name0, "<f8", " ",
                 1.0, "gbxindex", ngbxs),
        mom1zarr(store, maxchunk, name1, "<f8", "g",
                 dlc::MASS0grams, "gbxindex", ngbxs),
        mom2zarr(store, maxchunk, name2, "<f8", "g^2",
                 pow(dlc::MASS0grams, 2.0), "gbxindex", ngbxs){};
};

struct MomentsStorages : MassMom012Storages
{
  MomentsStorages(FSStore &store, const unsigned int maxchunk,
              const unsigned int ngbxs)
      : MassMom012Storages(store, maxchunk, ngbxs,
                    "mom0", "mom1", "mom2"){};
};

struct RainMomentsStorages : MassMom012Storages
{
  RainMomentsStorages(FSStore &store, const unsigned int maxchunk,
                  const unsigned int ngbxs)
      : MassMom012Storages(store, maxchunk, ngbxs,
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

#endif // MASSMOMENTSSTORAGE_HPP