/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: gbx_bounds_from_binary.hpp
 * Project: initialise
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions for reading gridbox boundaries from
 * a binary file (used to then create a
 * map from gbxindexes to gridbox boundaries
 * for CLEO SDM e.g. a CartesianMaps)
 */

#ifndef LIBS_INITIALISE_GBX_BOUNDS_FROM_BINARY_HPP_
#define LIBS_INITIALISE_GBX_BOUNDS_FROM_BINARY_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "initialise/readbinary.hpp"

/* holds vectors containing gridbox indexes and their
corresponding [zmin, zmax, zmin, xmax, ymin, ymax]
coordinate boundaries which are read from gridfile
and used in construction of GridboxMaps */
struct GbxBoundsFromBinary {
 private:
  void is_ngbxs_compatible(const size_t ngbxs) const;
  void is_nspacedims_compatible(const unsigned int nspacedims) const;

  bool check_0Dmodel_gbxbounds() const;
  bool check_1Dmodel_gbxbounds() const;
  bool check_2Dmodel_gbxbounds() const;
  bool check_3Dmodel_gbxbounds() const;

  size_t find_idx_in_gbxidxs(const unsigned int idx) const;

 public:
  std::vector<size_t> ndims;          // number of gridboxes in [coord3, coord1, coord2] dimensions
  std::vector<unsigned int> gbxidxs;  // gridbox indexes
  std::vector<double> gbxbounds;      // corresponding [coord3 {l, u}, coord1 {l, u}, coord2 {l, u}]
                                      // lower and upper coordinate boundaries

  GbxBoundsFromBinary(const size_t ngbxs, const unsigned int nspacedims,
                      const std::filesystem::path grid_filename);

  /* returns coord3 {lower, upper} gridbox bounds
  from position in gbxbounds vector which corresponds
  to position in gbxidxs where gbxidx = idx */
  Kokkos::pair<double, double> get_coord3gbxbounds(const unsigned int idx) const;

  /* returns coord1 {lower, upper} gridbox bounds
  from position in gbxbounds vector which corresponds
  to position in gbxidxs where gbxidx = idx */
  Kokkos::pair<double, double> get_coord1gbxbounds(const unsigned int idx) const;

  /* returns coord2 {lower, upper} gridbox bounds
  from position in gbxbounds vector which corresponds
  to position in gbxidxs where gbxidx = idx */
  Kokkos::pair<double, double> get_coord2gbxbounds(const unsigned int idx) const;

  /* calculates horizontal (x-y planar) area of gridbox
  using boundaries corresponding to gridbox with gbxidx=idx. */
  double gbxarea(const unsigned int idx) const;

  /* calculates volume of gridbox using boundaries
  corresponding to gridbox with gbxidx=idx. */
  double gbxvol(const unsigned int idx) const;

  /* returns total number of gridboxes = product of dimentions */
  size_t get_ngbxs() const { return ndims.at(0) * ndims.at(1) * ndims.at(2); }
};

#endif  // LIBS_INITIALISE_GBX_BOUNDS_FROM_BINARY_HPP_
