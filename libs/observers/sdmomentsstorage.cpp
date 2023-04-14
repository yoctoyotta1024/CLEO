// Author: Clara Bayley
// File: "sdmomentsstorage.cpp"
/* structs/classes to create an observer for the moments of
the superdroplet mass distribution that writes
into 1 dimensional array(s) 
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#include "sdmomentsstorage.hpp"

double mass0thmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx)
/* calculates the 0th moment of the (real) droplet mass distirbution
given by the superdrops in the span passed as an argument */
{
  double massmoment0 = 0.0;
  
  for (const auto &SDinGBx : span4SDsinGBx)
  {
    //massmoment0 += SDinGBx.superdrop.eps * SDinGBx.superdrop.mass();
    massmoment0 += SDinGBx.superdrop.eps;
  }
  
  return massmoment0;
}

double massnthmoment(const std::span<SuperdropWithGbxindex> span4SDsinGBx,
                      const double nth_moment)
/* calculates the nth moment of the (real) droplet mass distirbution
given by the superdrops in the span passed as an argument */
{
  double massmoment = 0.0;
  
  for (const auto &SDinGBx : span4SDsinGBx)
  {
    massmoment += SDinGBx.superdrop.eps * pow(SDinGBx.superdrop.mass(), nth_moment);
  }

  return massmoment;
}