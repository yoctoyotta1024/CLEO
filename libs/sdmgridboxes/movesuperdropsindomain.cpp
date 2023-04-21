
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
  if (at_domainboundary(index, 1, gbxmaps.ndims.at(0))) // at lower z edge of domain
  {
    superdrop.coord3 = coord3_beyondzdown(superdrop.coord3);
  }
  
  return gbxmaps.get_neighbour_zdown(index);
};

unsigned int zup(const Maps4GridBoxes &gbxmaps,
                 const unsigned int index,
                 Superdrop &superdrop)
/* function to update superdrop coord3 and return
sd_gbxindex of neighbouring gridbox in upwards z direction */
{
  if (at_domainboundary(index + 1, 1, gbxmaps.ndims.at(0))) // at upper z edge of domain
  {
    superdrop.coord3 = coord3_beyondzup(superdrop.coord3);
  }

  return gbxmaps.get_neighbour_zup(index);
};

unsigned int xbehind(const Maps4GridBoxes &gbxmaps,
                     const unsigned int index,
                     Superdrop &superdrop)
{
  const unsigned int increment = gbxmaps.ndims.at(0);
  if (at_domainboundary(index, increment, gbxmaps.ndims.at(1))) // at lower x edge of domain
  {
    superdrop.coord1 = coord1_beyondxbehind(superdrop.coord1);
  }

  return gbxmaps.get_neighbour_xbehind(index);
};

unsigned int xinfront(const Maps4GridBoxes &gbxmaps,
                      const unsigned int index,
                      Superdrop &superdrop)
{
  const unsigned int increment = gbxmaps.ndims.at(0);
  if (at_domainboundary(index + increment, increment, gbxmaps.ndims.at(1))) // at upper x edge of domain
  {
    superdrop.coord1 = coord1_beyondxinfront(superdrop.coord1);
  }

  return gbxmaps.get_neighbour_xinfront(index);
};

unsigned int yleft(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop)
{
  const unsigned int increment = gbxmaps.ndims.at(0) * gbxmaps.ndims.at(1); // no. gridboxes in z direction * no. gridboxes in x direction
  if (at_domainboundary(index, increment, gbxmaps.ndims.at(2))) // at lower y edge of domain
  {
    superdrop.coord2 = coord2_beyondyleft(superdrop.coord2);
  }

  return gbxmaps.get_neighbour_yleft(index);
};

unsigned int yright(const Maps4GridBoxes &gbxmaps,
                    const unsigned int index,
                    Superdrop &superdrop)
{
  const unsigned int increment = gbxmaps.ndims.at(0) * gbxmaps.ndims.at(1); // no. gridboxes in z direction * no. gridboxes in x direction
  if (at_domainboundary(index + increment, increment, gbxmaps.ndims.at(2))) // at upper y edge of domain
  {
    superdrop.coord2 = coord2_beyondyright(superdrop.coord2);
  }

  return gbxmaps.get_neighbour_yright(index);
};