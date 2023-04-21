
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
  if superdrop.coord3 <= domainlowerbound
  {
    superdrop.coord3 = coord3_beyondzdown();
  }
  
  return gbxmaps.get_neighbour_zdown(index);
};

unsigned int zup(const Maps4GridBoxes &gbxmaps,
                 const unsigned int index,
                 Superdrop &superdrop)
/* function to update superdrop coord3 and return
sd_gbxindex of neighbouring gridbox in upwards z direction */
{
  if superdrop.coord3 > domainupperbound 
  {
    superdrop.coord3 = coord3_beyondzup();
  }

  return gbxmaps.get_neighbour_zup(index);
};

unsigned int xbehind(const Maps4GridBoxes &gbxmaps,
                     const unsigned int index,
                     Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_xbehind(index);
};

unsigned int xinfront(const Maps4GridBoxes &gbxmaps,
                      const unsigned int index,
                      Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_xinfront(index);
};

unsigned int yleft(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_yleft(index);
};

unsigned int yright(const Maps4GridBoxes &gbxmaps,
                    const unsigned int index,
                    Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_yright(index);
};