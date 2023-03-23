// Author: Clara Bayley
// File: "observer_moments.hpp"
/* structs/classes to create an observer for the moments of
the superdroplet mass distribution that writes
into 1 dimensional array(s) 
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#ifndef OBSERVER_MOMENTS_HPP 
#define OBSERVER_MOMENTS_HPP 

#include <span>

#include "superdrop_solver/superdrop.hpp"
#include "./observer_singlevariable.hpp"

struct SDMomentsStorage
{
  TwoDStorage<double> massmoment0zarr;

  SDMomentsStorage(FSStore &store, const unsigned int maxcsize,
                   const unsigned int ngridboxes)
  : massmoment0zarr(store, maxcsize, "massmoment0", "<f8", " ", 1, ngridboxes)
  {};

};

double mass0thmoment(const std::span<SuperdropWithGridbox> span4SDsinGBx);
/* calculates the 0th moment of the (real) droplet mass distirbution
given by the superdrops in the span passed as an argument */

#endif // OBSERVER_MOMENTS_HPP 