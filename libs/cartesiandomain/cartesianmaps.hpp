/*
 * ----- CLEO -----
 * File: cartesianmaps.hpp
 * Project: cartesiandomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions related to creating and using maps to convert
 * between a gridbox indexes and domain coordinates for a
 * cartesian C grid
 */

#ifndef CARTESIANMAPS_HPP
#define CARTESIANMAPS_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"

// TODO: host versions? (maps are in exec mem space)

namespace dlc = dimless_constants;

struct CartesianMaps
/* type satisfying GridboxMaps concept specifically for gridboxes
defined on in a cartesian C grid with equal area and volume for each 
gridbox. coord[X]bounds (for X = 1, 2, 3, corresponding to x, y, z) 
map between gbxindexes and gridbox boundaries in a cartiesian domain.
The keys of each map are the gridbox indexes. The corresponding value
of each bounds map is that gridbox's {lower boundary, upper boundary}.
to_[direction]_coord[X]nghbr (for direction = back, forward)
map from a given gbxidx to the gbxidx of a neighbouring gridbox
in that direction */
{
private:
  /* additional gridbox / domain information */
  viewd_ndims ndims; // dimensions (ie. no. gridboxes) in [coord3, coord1, coord2] directions
  double gbxareas;   // horizontal (x-y planar) area of all gridboxes
  double gbxvolumes; // volume of all gridboxes

public:
  /* maps from gbxidx to {lower, upper} coords of gridbox boundaries */
  kokkos_pairmap to_coord3bounds;
  kokkos_pairmap to_coord1bounds;
  kokkos_pairmap to_coord2bounds;

  /* maps from gbxidx to gbxindx of front / back neighbour */
  kokkos_uintmap to_back_coord3nghbr;
  kokkos_uintmap to_forward_coord3nghbr;
  kokkos_uintmap to_back_coord1nghbr;
  kokkos_uintmap to_forward_coord1nghbr;
  kokkos_uintmap to_back_coord2nghbr;
  kokkos_uintmap to_forward_coord2nghbr;

  CartesianMaps(const size_t ngbxs)
  {

  }

  KOKKOS_INLINE_FUNCTION CartesianMaps() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~CartesianMaps() = default; // Kokkos requirement for a (dual)View

  void set_ndims(viewd_ndims d_ndims)
  {
    ndims = d_ndims;
  }

  void set_gbxarea(const double iarea)
  {
    gbxareas = iarea;
  }

  void set_gbxvolume(const double ivolume)
  {
    gbxvolumes = ivolume;
  }

  KOKKOS_INLINE_FUNCTION
  double get_gbxarea(const unsigned int gbxidx) const
  {
    return gbxareas;
  }

  KOKKOS_INLINE_FUNCTION
  double get_gbxvolume(const unsigned int gbxidx) const
  {
    return gbxvolumes;
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  coord3bounds(const unsigned int gbxidx) const
  /* returns {lower bound, upper bound}  in coord3
  (z) direction of gridbox with index 'gbxidx'
  on device */
  {
    const auto i(to_coord3bounds.find(gbxidx)); // index in map of key 'gbxidx'

    return to_coord3bounds.value_at(i); // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  coord1bounds(const unsigned int gbxidx) const
  /* returns {lower bound, upper bound}  in coord1
  (x) direction of gridbox with index 'gbxidx'
  on device */
  {
    const auto i(to_coord1bounds.find(gbxidx)); // index in map of key 'gbxidx'

    return to_coord1bounds.value_at(i); // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  coord2bounds(const unsigned int gbxidx) const
  /* returns {lower bound, upper bound}  in coord2
  (y) direction of gridbox with index 'gbxidx'
  on device */
  {
    const auto i(to_coord2bounds.find(gbxidx)); // index in map of key 'gbxidx'

    return to_coord2bounds.value_at(i); // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int coord3backward(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the backwards coord3, ie. downwards z, direction */
  {
    const auto i(to_back_coord3nghbr.find(gbxindex)); // index in map of key 'gbxidx'

    return to_back_coord3nghbr.value_at(i); // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int coord3forward(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards coord3, ie. upwards z, direction */
  {
    const auto i(to_forward_coord3nghbr.find(gbxindex)); // index in map of key 'gbxidx'

    return to_forward_coord3nghbr.value_at(i); // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int coord1backward(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the backwards coord1, ie. into page x, direction */
  {
    const auto i(to_back_coord1nghbr.find(gbxindex)); // index in map of key 'gbxidx'

    return to_back_coord1nghbr.value_at(i); // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int coord1forward(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards coord1, ie. out of page x, direction */
  {
    const auto i(to_forward_coord1nghbr.find(gbxindex)); // index in map of key 'gbxidx'

    return to_forward_coord1nghbr.value_at(i); // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int coord2backward(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the backwards coord2, ie. left y, direction */
  {
    const auto i(to_back_coord2nghbr.find(gbxindex)); // index in map of key 'gbxidx'

    return to_back_coord2nghbr.value_at(i); // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int coord2forward(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards coord2, ie. right y, direction */
  {
    const auto i(to_forward_coord2nghbr.find(gbxindex)); // index in map of key 'gbxidx'

    return to_forward_coord2nghbr.value_at(i); // value returned by map at index i
  }
};

#endif // CARTESIANMAPS_HPP