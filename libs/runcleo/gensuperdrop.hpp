/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: gensuperdrop.hpp
 * Project: runcleo
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * classes and templated functions required by RunCLEO to create a view of superdroplets
 * (on device) using some initial conditions
 */

#ifndef LIBS_RUNCLEO_GENSUPERDROP_HPP_
#define LIBS_RUNCLEO_GENSUPERDROP_HPP_

#include <Kokkos_Core.hpp>
#include <array>
#include <memory>

#include "../kokkosaliases.hpp"
#include "initialise/initialconditions.hpp"
#include "superdrops/superdrop.hpp"

/**
 * @brief Struct that holds data for the initial conditions of super-droplets.
 *
 * This struct holds data for the initial conditions of some properties of super-droplets
 * and provides an operator which returns a super-droplet generated using them.
 */
class GenSuperdrop {
 private:
  size_t maxnsupers; /**< Maximum allowed superdroplets in the domain (total extent of view). */
  unsigned int nspacedims; /**< Number of spatial dimensions. */
  InitSupersData initdata; /**< instance of InitSupersData for initialising superdrops. */

  /**
   * @brief Returns initial spatial coordinates of the kk'th super-droplet.
   *
   * A coordinate is only copied from the corresponding coords vector if that coordinate is
   * consistent with the number of spatial dimensions of the model. Otherwise, the coordinate is set
   * to 0.0. For example, if the model is 1-D, only coord3 is obtained from the initial data vector;
   * coord1 and coord2 are set to 0.0.
   *
   * @param kk The index of the super-droplet in the initial
   * data (0 <= kk < total number of super-droplets).
   * @return An array containing the kk'th superdroplet's spatial
   * coordinates (coord3, coord1, coord2).
   */
  std::array<double, 3> coords_at(const unsigned int kk) const;

  /* helper function to return a superdroplet's attributes
  at position kk in the initial conditions data. All
  superdroplets created with same solute properties */
  /**
   * @brief Function returns a superdroplet's attributes
   * from position 'kk' in the initial conditions data. All
   * super-droplets have the same solute properties.
   *
   * @param kk The index of the super-droplet in the initial
   * data (0 <= kk < total number of super-droplets).
   * @return The attributes of the superdrop from index 'kk'.
   */
  SuperdropAttrs attrs_at(const unsigned int kk) const;

 public:
  /**
   * @brief Constructs a GenSuperdrop instance.
   *
   * This constructor initializes a GenSuperdrop instance using the provided
   * SuperdropInitConds instance to fetch initial data.
   *
   * @tparam SuperdropInitConds The type of the data for super-droplets' initial conditions.
   * @param sdic The instance of the data for super-droplets' initial conditions.
   */
  template <typename SuperdropInitConds>
  explicit GenSuperdrop(const SuperdropInitConds &sdic)
      : maxnsupers(sdic.get_maxnsupers()),
        nspacedims(sdic.get_nspacedims()),
        initdata(sdic.fetch_data()) {}

  auto get_maxnsupers() { return maxnsupers; }

  /**
   * @brief Generate a super-droplet using initial data for the kk'th superdrop.
   *
   * This function returns a superdrop generated from the specified position
   * 'kk' in the initial conditions data.
   *
   * @param kk The index of the superdrop to generate.
   * @return The generated super-droplet.
   */
  Superdrop operator()(const unsigned int kk) const;
};

#endif  // LIBS_RUNCLEO_GENSUPERDROP_HPP_
