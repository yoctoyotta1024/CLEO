
// Author: Clara Bayley
// File: "movesuperdropsindomain.cpp"
/* Implementation of functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#include "./movesuperdropsindomain.hpp"

unsigned int zdown(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop)
/* function to update superdrop coord3 and return
sd_gbxindex of neighbouring gridbox in downards z direction */
{
  const unsigned int nghbour(gbxmaps.get_neighbour_zdown(index));

  if (at_domainboundary(index, 1, gbxmaps.ndims.at(0))) // at lower z edge of domain
  {
    const double lim1 = gbxmaps.get_bounds_z(nghbour).second; // upper lim of backward nghbour
    const double lim2 = gbxmaps.get_bounds_z(index).first; // lower lim of gbx
    superdrop.coord3 = coord3_beyondz(superdrop.coord3, lim1, lim2);
  }
  
  return nghbour; // gbxindex of zdown_neighbour
};

unsigned int zup(const Maps4GridBoxes &gbxmaps,
                 const unsigned int index,
                 Superdrop &superdrop)
/* function to update superdrop coord3 and return
sd_gbxindex of neighbouring gridbox in upwards z direction */
{
  const unsigned int nghbour(gbxmaps.get_neighbour_zup(index));

  if (at_domainboundary(index + 1, 1, gbxmaps.ndims.at(0))) // at upper z edge of domain
  {
    const double lim1 = gbxmaps.get_bounds_z(nghbour).first; // lower lim of forward nghbour
    const double lim2 = gbxmaps.get_bounds_z(index).second; // upper lim of gbx
    superdrop.coord3 = coord3_beyondz(superdrop.coord3, lim1, lim2);
  }

  return nghbour; // gbxindex of zup_neighbour
};

unsigned int xbehind(const Maps4GridBoxes &gbxmaps,
                     const unsigned int index,
                     Superdrop &superdrop)
{
  const unsigned int nghbour(gbxmaps.get_neighbour_xbehind(index));
  
  const unsigned int increment = gbxmaps.ndims.at(0);
  if (at_domainboundary(index, increment, gbxmaps.ndims.at(1))) // at lower x edge of domain
  {
    const double lim1 = gbxmaps.get_bounds_x(nghbour).second; // upper lim of backward nghbour
    const double lim2 = gbxmaps.get_bounds_x(index).first; // lower lim of gbx
    superdrop.coord1 = coord1_beyondx(superdrop.coord1, lim1, lim2);
  }

  return nghbour; // gbxindex of xbehind_neighbour
};

unsigned int xinfront(const Maps4GridBoxes &gbxmaps,
                      const unsigned int index,
                      Superdrop &superdrop)
{
  const unsigned int nghbour(gbxmaps.get_neighbour_xinfront(index));

  const unsigned int increment = gbxmaps.ndims.at(0);
  if (at_domainboundary(index + increment, increment, gbxmaps.ndims.at(1))) // at upper x edge of domain
  {
    const double lim1 = gbxmaps.get_bounds_x(nghbour).first; // lower lim of forward nghbour
    const double lim2 = gbxmaps.get_bounds_x(index).second; // upper lim of gbx
    superdrop.coord1 = coord1_beyondx(superdrop.coord1, lim1, lim2);
  }

  return nghbour;
};

unsigned int yleft(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop)
{
  const unsigned int nghbour(gbxmaps.get_neighbour_yleft(index));

  const unsigned int increment = gbxmaps.ndims.at(0) * gbxmaps.ndims.at(1); // no. gridboxes in z direction * no. gridboxes in x direction
  if (at_domainboundary(index, increment, gbxmaps.ndims.at(2))) // at lower y edge of domain
  {
    const double lim1 = gbxmaps.get_bounds_y(nghbour).second; // upper lim of backward nghbour
    const double lim2 = gbxmaps.get_bounds_y(index).first; // lower lim of gbx
    superdrop.coord2 = coord2_beyondy(superdrop.coord2, lim1, lim2);
  }

  return nghbour // gbxindex of yleft_neighbour
};

unsigned int yright(const Maps4GridBoxes &gbxmaps,
                    const unsigned int index,
                    Superdrop &superdrop)
{
  const unsigned int nghbour(gbxmaps.get_neighbour_yright(index));

  const unsigned int increment = gbxmaps.ndims.at(0) * gbxmaps.ndims.at(1); // no. gridboxes in z direction * no. gridboxes in x direction
  if (at_domainboundary(index + increment, increment, gbxmaps.ndims.at(2))) // at upper y edge of domain
  {
    const double lim1 = gbxmaps.get_bounds_y(nghbour).first; // lower lim of forward nghbour
    const double lim2 = gbxmaps.get_bounds_y(index).second; // upper lim of gbx
    superdrop.coord2 = coord2_beyondy(superdrop.coord2, lim1, lim2);
  }

  return nghbour;
};