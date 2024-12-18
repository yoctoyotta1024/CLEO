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
 * Last Modified: Friday 3rd May 2024
 * Modified By: CB
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
// TODO(CB): use domain_decomposition instead of ndims
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
  viewd_ndims ndims;  // dimensions (ie. no. gridboxes) in [coord3, coord1, coord2] directions //
                      // TODO(CB): delete

 public:
  /* initialise maps with hint for their capacity
  (ie. total number of keys). Leave values for maps,
  for ndims, gbxareas and gbxvols undefined
  upon construction (e.g. null ptr for ndims) */
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
        ndims("ndims") {}

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

  /* insert 1 value into to_areas map at key = idx with value=area */
  void insert_gbxarea(const unsigned int idx, double area) const { to_areas.insert(idx, area); }

  /* insert 1 value into to_volumes map at key = idx with value=volume */
  void insert_gbxvolume(const unsigned int idx, double volume) const {
    to_volumes.insert(idx, volume);
  }

  /* insert 1 value into to_coord3bounds
  map at key = idx with value=bounds TODO(CB): WIP delete function (and analogous ones) */
  void insert_coord3bounds(const unsigned int idx, Kokkos::pair<double, double> bounds) {
    /* parallel for for 1 value so that execution of
    insert occurs on device if necessary */
    Kokkos::parallel_for(
        "cb3s", 1,
        KOKKOS_CLASS_LAMBDA(const unsigned int i) { to_coord3bounds.insert(idx, bounds); });
  }

  /* insert 1 value into to_coord1bounds
  map at key = idx with value=bounds */
  void insert_coord1bounds(const unsigned int idx, Kokkos::pair<double, double> bounds) {
    /* parallel for for 1 value so that execution of
    insert occurs on device if necessary */
    Kokkos::parallel_for(
        "cb1s", 1,
        KOKKOS_CLASS_LAMBDA(const unsigned int i) { to_coord1bounds.insert(idx, bounds); });
  }

  /* insert 1 value into to_coord2bounds
  map at key = idx with value=bounds */
  void insert_coord2bounds(const unsigned int idx, Kokkos::pair<double, double> bounds) {
    /* parallel for for 1 value so that execution of
    insert occurs on device if necessary */
    Kokkos::parallel_for(
        "cb2s", 1,
        KOKKOS_CLASS_LAMBDA(const unsigned int i) { to_coord2bounds.insert(idx, bounds); });
  }

  /* insert indexes of neightbouring gridboxes into
  back and forward neighbour maps for key = idx given
  neighbours in pair {back, forward},
  i.e. back=nghbrs.first and forward=nghbrs.second */
  void insert_coord3nghbrs(const unsigned int idx,
                           Kokkos::pair<unsigned int, unsigned int> nghbrs) {
    /* parallel for for 1 value so that execution of
    insert occurs on device if necessary */
    Kokkos::parallel_for(
        "nghbr3s", 1, KOKKOS_CLASS_LAMBDA(const unsigned int i) {
          to_back_coord3nghbr.insert(idx, nghbrs.first);
          to_forward_coord3nghbr.insert(idx, nghbrs.second);
        });
  }

  /* insert indexes of neightbouring gridboxes into
  back and forward neighbour maps for key = idx given
  neighbours in pair {back, forward},
  i.e. back=nghbrs.first and forward=nghbrs.second */
  void insert_coord1nghbrs(const unsigned int idx,
                           Kokkos::pair<unsigned int, unsigned int> nghbrs) {
    /* parallel for for 1 value so that execution of
    insert occurs on device if necessary */
    Kokkos::parallel_for(
        "nghbr1s", 1, KOKKOS_CLASS_LAMBDA(const unsigned int i) {
          to_back_coord1nghbr.insert(idx, nghbrs.first);
          to_forward_coord1nghbr.insert(idx, nghbrs.second);
        });
  }

  /* insert indexes of neightbouring gridboxes into
  back and forward neighbour maps for key = idx given
  neighbours in pair {back, forward},
  i.e. back=nghbrs.first and forward=nghbrs.second */
  void insert_coord2nghbrs(const unsigned int idx,
                           Kokkos::pair<unsigned int, unsigned int> nghbrs) {
    /* parallel for for 1 value so that execution of
    insert occurs on device if necessary */
    Kokkos::parallel_for(
        "nghbr2s", 1, KOKKOS_CLASS_LAMBDA(const unsigned int i) {
          to_back_coord2nghbr.insert(idx, nghbrs.first);
          to_forward_coord2nghbr.insert(idx, nghbrs.second);
        });
  }

  /* copies of h_ndims to ndims, possibly into device memory */
  void set_ndims_via_copy(const viewd_ndims::HostMirror h_ndims) {
    Kokkos::deep_copy(ndims, h_ndims);
  }

  /* returns model dimensions ie. number of gridboxes
  along [coord3, coord1, coord2] directions for use on
  host. deep copy is made if gbxmaps ndims is on device */
  viewd_ndims::HostMirror get_ndims_hostcopy() const {
    auto h_ndims = Kokkos::create_mirror_view(ndims);  // mirror in case ndims is on device memory
    Kokkos::deep_copy(h_ndims, ndims);

    return h_ndims;
  }

  /* returns model dimensions ie. number of gridboxes
  along [coord3, coord1, coord2] directions */
  KOKKOS_INLINE_FUNCTION
  viewd_ndims get_ndims() const { return ndims; }

  /* returns model dimensions ie. number of
  gridboxes along d'th direction, where:
  ndims(d=0) = coord3
  ndims(d=1) = coord1
  ndims(d=2) = coord2 */
  KOKKOS_INLINE_FUNCTION
  size_t get_ndim(const unsigned int d) const { return ndims(d); }

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

  void create_decomposition(std::vector<size_t> ndims, double gridbox_z_size, double gridbox_x_size,
                            double gridbox_y_size) {
    domain_decomposition.create(ndims, gridbox_z_size, gridbox_x_size, gridbox_y_size);
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
  unsigned int get_local_bounding_gridbox(const unsigned int gbxindex, double& coord3,
                                          double& coord1, double& coord2) const;
};

#endif  // LIBS_CARTESIANDOMAIN_CARTESIANMAPS_HPP_
