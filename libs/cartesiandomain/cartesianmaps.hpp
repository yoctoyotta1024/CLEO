/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesianmaps.hpp
 * Project: cartesiandomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions related to creating and using maps to convert
 * between a gridbox indexes and domain coordinates for a
 * cartesian C grid
 */

#ifndef LIBS_CARTESIANDOMAIN_CARTESIANMAPS_HPP_
#define LIBS_CARTESIANDOMAIN_CARTESIANMAPS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include <array>
#include <stdexcept>
#include <vector>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesian_decomposition.hpp"
#include "cartesiandomain/doubly_periodic_domain.hpp"
#include "initialise/gbx_bounds_from_binary.hpp"

namespace dlc = dimless_constants;

/* type satisfying GridboxMaps concept specifically for gridboxes
defined on in a cartesian C grid with equal area and volume for each
gridbox. coord[X]bounds (for X = 1, 2, 3, corresponding to x, y, z)
map between gbxindexes and gridbox boundaries in a cartiesian domain.
The keys of each map are the gridbox indexes. The corresponding value
of each bounds map is that gridbox's {lower boundary, upper boundary}.
to_[direction]_coord[X]nghbr (for direction = back, forward)
map from a given gbxidx to the gbxidx of a neighbouring gridbox
in that direction */
// TODO(CB): use domain_decomposition instead of global_ndims
struct CartesianMaps {
 private:
  CartesianDecomposition domain_decomposition;
  bool is_decomp;

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

  /* additional gridbox / domain information */
  kokkos_dblmaph to_areas;  // map from gbxindex to horizontal (x-y planar) area of gridbox on host
  kokkos_dblmaph to_volumes;  // map from gbxindex to volume of gridbox on host
  viewd_ndims global_ndims;   // entire domain no. gridboxes in [coord3, coord1, coord2] dimensions

 public:
  /* initialise maps without capacity. Note this leave values for maps,
  for e.g. for global_ndims, gbxareas and gbxvols undefined upon construction */
  explicit CartesianMaps()
      : is_decomp(false),
        to_coord3bounds(kokkos_pairmap(0)),
        to_coord1bounds(kokkos_pairmap(0)),
        to_coord2bounds(kokkos_pairmap(0)),
        to_back_coord3nghbr(kokkos_uintmap(0)),
        to_forward_coord3nghbr(kokkos_uintmap(0)),
        to_back_coord1nghbr(kokkos_uintmap(0)),
        to_forward_coord1nghbr(kokkos_uintmap(0)),
        to_back_coord2nghbr(kokkos_uintmap(0)),
        to_forward_coord2nghbr(kokkos_uintmap(0)),
        to_areas(kokkos_dblmaph(0)),
        to_volumes(kokkos_dblmaph(0)),
        global_ndims("global_ndims") {}

  /* copy host version of to_coord3bounds to gridbox maps (possibly in device memory) */
  void set_coord3bounds_via_copy(const kokkos_pairmap::HostMirror h_to_coord3bounds) {
    to_coord3bounds.create_copy_view(h_to_coord3bounds);
  }

  /* copy host version of to_coord1bounds to gridbox maps (possibly in device memory) */
  void set_coord1bounds_via_copy(const kokkos_pairmap::HostMirror h_to_coord1bounds) {
    to_coord1bounds.create_copy_view(h_to_coord1bounds);
  }

  /* copy host version of to_coord2bounds to gridbox maps (possibly in device memory) */
  void set_coord2bounds_via_copy(const kokkos_pairmap::HostMirror h_to_coord2bounds) {
    to_coord2bounds.create_copy_view(h_to_coord2bounds);
  }

  /* copy host version of to_back_coord3nghbr to gridbox maps (possibly in device memory) */
  void set_back_coord3nghbr_via_copy(const kokkos_uintmap::HostMirror h_to_back_coord3nghbr) {
    to_back_coord3nghbr.create_copy_view(h_to_back_coord3nghbr);
  }

  /* copy host version of to_forward_coord3nghbr to gridbox maps (possibly in device memory) */
  void set_forward_coord3nghbr_via_copy(const kokkos_uintmap::HostMirror h_to_forward_coord3nghbr) {
    to_forward_coord3nghbr.create_copy_view(h_to_forward_coord3nghbr);
  }

  /* copy host version of to_back_coord1nghbr to gridbox maps (possibly in device memory) */
  void set_back_coord1nghbr_via_copy(const kokkos_uintmap::HostMirror h_to_back_coord1nghbr) {
    to_back_coord1nghbr.create_copy_view(h_to_back_coord1nghbr);
  }

  /* copy host version of to_forward_coord1nghbr to gridbox maps (possibly in device memory) */
  void set_forward_coord1nghbr_via_copy(const kokkos_uintmap::HostMirror h_to_forward_coord1nghbr) {
    to_forward_coord1nghbr.create_copy_view(h_to_forward_coord1nghbr);
  }

  /* copy host version of to_back_coord2nghbr to gridbox maps (possibly in device memory) */
  void set_back_coord2nghbr_via_copy(const kokkos_uintmap::HostMirror h_to_back_coord2nghbr) {
    to_back_coord2nghbr.create_copy_view(h_to_back_coord2nghbr);
  }

  /* copy host version of to_forward_coord2nghbr to gridbox maps (possibly in device memory) */
  void set_forward_coord2nghbr_via_copy(const kokkos_uintmap::HostMirror h_to_forward_coord2nghbr) {
    to_forward_coord2nghbr.create_copy_view(h_to_forward_coord2nghbr);
  }

  void set_gbxareas_map(const kokkos_dblmaph i_to_areas) { to_areas = i_to_areas; }

  void set_gbxvolumes_map(const kokkos_dblmaph i_to_volumes) { to_volumes = i_to_volumes; }

  /* copies of h_global_ndims to global_ndims, possibly into device memory */
  void set_global_ndims_via_copy(const viewd_ndims::HostMirror h_global_ndims) {
    Kokkos::deep_copy(global_ndims, h_global_ndims);
  }

  /* returns model dimensions ie. number of gridboxes along [coord3, coord1, coord2]
  directions for use on host. deep copy is made if gbxmaps global_ndims is in device memory */
  viewd_ndims::HostMirror get_global_ndims_hostcopy() const {
    auto h_global_ndims = Kokkos::create_mirror_view(global_ndims);
    Kokkos::deep_copy(h_global_ndims, global_ndims);

    return h_global_ndims;
  }

  /* returns model dimensions ie. number of gridboxes
  along [coord3, coord1, coord2] directions */
  KOKKOS_INLINE_FUNCTION
  viewd_ndims get_global_ndims() const { return global_ndims; }

  /* returns model dimensions ie. number of
  gridboxes along d'th direction, where:
  global_ndims(d=0) = coord3
  global_ndims(d=1) = coord1
  global_ndims(d=2) = coord2 */
  KOKKOS_INLINE_FUNCTION
  size_t get_global_ndim(const unsigned int d) const { return global_ndims(d); }

  /* on host device, throws error if maps are not all
  the same size, else returns size of maps */
  size_t maps_size() const;

  /* returns volume of gridbox with index 'gbxidx' on host */
  double get_gbxvolume(const unsigned int gbxidx) const {
    const auto i = to_volumes.find(gbxidx);  // index in map of key 'gbxindex'

    return to_volumes.value_at(i);  // value returned by map at index i
  }

  /* returns horizontal (x-y planar) area of gridbox with index 'gbxidx' on host */
  double get_gbxarea(const unsigned int gbxidx) const {
    const auto i = to_areas.find(gbxidx);  // index in map of key 'gbxindex'

    return to_areas.value_at(i);  // value returned by map at index i
  }

  /* returns {lower bound, upper bound}  in coord3
  (z) direction of gridbox with index 'gbxidx'
  on device */
  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double> coord3bounds(const unsigned int gbxidx) const {
    const auto i = to_coord3bounds.find(gbxidx);  // index in map of key 'gbxidx'

    return to_coord3bounds.value_at(i);  // value returned by map at index i
  }

  /* returns {lower bound, upper bound}  in coord1
  (x) direction of gridbox with index 'gbxidx'
  on device */
  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double> coord1bounds(const unsigned int gbxidx) const {
    const auto i = to_coord1bounds.find(gbxidx);  // index in map of key 'gbxidx'

    return to_coord1bounds.value_at(i);  // value returned by map at index i
  }

  /* returns {lower bound, upper bound}  in coord2
  (y) direction of gridbox with index 'gbxidx'
  on device */
  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double> coord2bounds(const unsigned int gbxidx) const {
    const auto i = to_coord2bounds.find(gbxidx);  // index in map of key 'gbxidx'

    return to_coord2bounds.value_at(i);  // value returned by map at index i
  }

  /* given gridbox index, return index of neighbouring
  gridbox in the backwards coord3, ie. downwards z, direction */
  KOKKOS_INLINE_FUNCTION
  unsigned int coord3backward(unsigned int gbxindex) const {
    const auto i = to_back_coord3nghbr.find(gbxindex);  // index in map of key 'gbxidx'

    return to_back_coord3nghbr.value_at(i);  // value returned by map at index i
  }

  /* given gridbox index, return index of neighbouring
  gridbox in the forwards coord3, ie. upwards z, direction */
  KOKKOS_INLINE_FUNCTION
  unsigned int coord3forward(unsigned int gbxindex) const {
    const auto i = to_forward_coord3nghbr.find(gbxindex);  // index in map of key 'gbxidx'

    return to_forward_coord3nghbr.value_at(i);  // value returned by map at index i
  }

  /* given gridbox index, return index of neighbouring
  gridbox in the backwards coord1, ie. into page x, direction */
  KOKKOS_INLINE_FUNCTION
  unsigned int coord1backward(unsigned int gbxindex) const {
    const auto i = to_back_coord1nghbr.find(gbxindex);  // index in map of key 'gbxidx'

    return to_back_coord1nghbr.value_at(i);  // value returned by map at index i
  }

  /* given gridbox index, return index of neighbouring
  gridbox in the forwards coord1, ie. out of page x, direction */
  KOKKOS_INLINE_FUNCTION
  unsigned int coord1forward(unsigned int gbxindex) const {
    const auto i = to_forward_coord1nghbr.find(gbxindex);  // index in map of key 'gbxidx'

    return to_forward_coord1nghbr.value_at(i);  // value returned by map at index i
  }

  /* given gridbox index, return index of neighbouring
  gridbox in the backwards coord2, ie. left y, direction */
  KOKKOS_INLINE_FUNCTION
  unsigned int coord2backward(unsigned int gbxindex) const {
    const auto i = to_back_coord2nghbr.find(gbxindex);  // index in map of key 'gbxidx'

    return to_back_coord2nghbr.value_at(i);  // value returned by map at index i
  }

  /* given gridbox index, return index of neighbouring
  gridbox in the forwards coord2, ie. right y, direction */
  KOKKOS_INLINE_FUNCTION
  unsigned int coord2forward(unsigned int gbxindex) const {
    const auto i = to_forward_coord2nghbr.find(gbxindex);  // index in map of key 'gbxidx'

    return to_forward_coord2nghbr.value_at(i);  // value returned by map at index i
  }

  void create_decomposition(std::vector<size_t> global_ndims, GbxBoundsFromBinary gfb ) {
    domain_decomposition.create(global_ndims, gfb);
    if (domain_decomposition.get_total_local_gridboxes() <
        domain_decomposition.get_total_global_gridboxes()) {
      is_decomp = true;
    }
  }

  const CartesianDecomposition& get_domain_decomposition() const { return domain_decomposition; }

  size_t get_total_global_ngridboxes() const {
    return domain_decomposition.get_total_global_gridboxes();
  }

  KOKKOS_FUNCTION
  size_t get_local_ngridboxes() const;

  // TODO(ALL): refactor once domain_decomposition.get_total_local_gridboxes() is a GPU function
  size_t get_local_ngridboxes_hostcopy() const {
    return domain_decomposition.get_total_local_gridboxes();
  }

  unsigned int global_to_local_gbxindex(size_t global_gridbox_index) const {
    return domain_decomposition.global_to_local_gridbox_index(global_gridbox_index);
  }

  KOKKOS_FUNCTION
  size_t local_to_global_gridbox_index(unsigned int local_gridbox_index, int process = -1) const;

  /* given coordinates, associated gxbindex is found. The coords may be updated too,
   * e.g. if the domain has a cyclic boundary condition and they therefore need to be corrected
   */
  KOKKOS_FUNCTION
  unsigned int get_local_bounding_gridbox_index(const unsigned int gbxindex, double& coord3,
                                          double& coord1, double& coord2) const;
};

#endif  // LIBS_CARTESIANDOMAIN_CARTESIANMAPS_HPP_
