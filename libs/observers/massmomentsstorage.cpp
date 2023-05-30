// Author: Clara Bayley
// File: "massmomentsstorage.cpp"
/* structs/classes to create an observer for the moments of
the superdroplet mass distribution that writes
into 1 dimensional array(s)
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#include "massmomentsstorage.hpp"

double massmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx,
                  const double nth_moment)
/* calculates the nth moment of the (real) droplet mass distirbution
given by the superdrops in the span passed as an argument */
{
  double nthmom = 0.0;
  for (const auto &SDinGBx : span4SDsinGBx)
  {
    nthmom += SDinGBx.superdrop.eps *
                  pow(SDinGBx.superdrop.mass(), nth_moment);
  }
  return nthmom;
}

double rainmassmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx,
                     const double nth_moment)
/* calculates the nth moment of the (real) raindroplet mass
distirbution given by all the superdrops which have radius >= rlim
in the span passed as an argument */
{
  const double rlim(40e-6/dlc::R0); // minimum dimless radius of a raindrop

  double nthmom = 0.0;
  for (const auto &SDinGBx : span4SDsinGBx)
  {
    if (SDinGBx.superdrop.radius >= rlim)
    {
      nthmom += SDinGBx.superdrop.eps *
                    pow(SDinGBx.superdrop.mass(), nth_moment);
    } 
  }
  return nthmom;
}