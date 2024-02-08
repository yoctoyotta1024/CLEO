/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: createsupers.hpp
 * Project: runcleo
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 8th February 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * file for structure(s) to create a view of
 * superdroplets (on device) using some
 * initial conditions
 */

#ifndef LIBS_RUNCLEO_CREATESUPERS_HPP_
#define LIBS_RUNCLEO_CREATESUPERS_HPP_

#include <array>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/sortsupers.hpp"
#include "initialise/initconds.hpp"
#include "superdrops/superdrop.hpp"

/**
 * @brief Struct that holds data for the initial conditions of super-droplets.
 *
 * This struct holds data for the initial conditions of some properties of super-droplets
 * and provides an operator which returns a super-droplet generated using them.
 */
class GenSuperdrop {
 private:
  unsigned int nspacedims;   ///< Number of spatial dimensions.
  std::shared_ptr<Superdrop::IDType::Gen> sdIdGen;   ///< Pointer to super-droplet ID generator.
  InitSupersData initdata;   ///< Data for initialising superdrops.

  /**
   * @brief Returns initial spatial coordinates of the kk'th super-droplet.
   *
   * A coordinate is only copied from the corresponding coords vector if that coordinate is
   * consistent with the number of spatial dimensions of the model. Otherwise, the coordinate is set
   * to 0.0. For example, if the model is 1-D, only coord3 is obtained from the initial data vector;
   * coord1 and coord2 are set to 0.0.
   *
   * @param kk The index of the superdroplet in the initial
   * data (0 <= kk < total number of superdrops).
   * @return An array containing the kk'th superdroplet's spatial
   * coordinates (coord3, coord1, coord2).
   */
  std::array<double, 3> coords_at(const unsigned int kk) const;

  /**
   * @brief Get the attributes of the kk'th superdrop.
   *
   * @param kk The index of the superdrop (0 <= kk < total number of superdrops).
   * @return The attributes of the superdrop at index 'kk'.
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
      : nspacedims(sdic.get_nspacedims()), sdIdGen(std::make_shared<Superdrop::IDType::Gen>()) {
    sdic.fetch_data(initdata);
  }

  /**
   * @brief Generate a superdrop using initial data for the kk'th superdrop.
   *
   * This function returns a superdrop generated using the initial conditions
   * for the kk'th superdrop.
   *
   * @param kk The index of the superdrop to generate.
   * @return The generated superdrop.
   */
  Superdrop operator()(const unsigned int kk) const;
};

/**
 * @brief Create a view of super-droplets in (device) memory.
 *
 * This function creates an ordered view of superdrops in device memory, where the number
 * of superdrops is specified by the parameter `totnsupers`. The superdrops are
 * ordered by the gridbox indexes and generated using a generator which uses
 * the initial conditions referenced by the `SuperdropInitConds` type.
 *
 * @tparam SuperdropInitConds The type of the super-droplets' initial conditions data.
 * @param sdic The instance of the super-droplets' initial conditions data.
 * @return A view of super-droplets in device memory.
 */
template <typename SuperdropInitConds>
viewd_supers create_supers(const SuperdropInitConds &sdic);

/**
 * @brief Return an initialised view of superdrops in device memory.
 *
 * This function initialises a view of superdrops in device memory by creating
 * a view on the device and copying a host mirror view that is initialised using
 * the `SuperdropInitConds` instance.
 *
 * @tparam SuperdropInitConds The type of the super-droplets' initial conditions data.
 * @param sdic The instance of the super-droplets' initial conditions data.
 * @return A view of superdrops in device memory.
 */
template <typename SuperdropInitConds>
viewd_supers initialise_supers(const SuperdropInitConds &sdic);

/**
 * @brief Return a mirror view of superdrops on host memory.
 *
 * This function initialises a mirror view of superdrops on host memory, using the super-droplet
 * generator instance `GenSuperdrop` to generate the kk'th super-droplet with their initial Gridbox
 * index, spatial coordinates, and attributes. Kokkos::parallel_for is used to perform parallel
 * initialisation of the mirror view, where each superdrop is generated by the provided `gen`
 * function object.
 *
 * The equivalent serial version of the Kokkos::parallel_for([...]) is:
 * @code
 * for (size_t kk(0); kk < totnsupers; ++kk)
 * {
 *  h_supers(kk) = gen(kk);
 * }
 * @endcode
 *
 * @param gen The generator function object used to initialize the superdrops.
 * @param supers The view of superdrops on device memory.
 * @return A mirror view of superdrops on host memory.
 */
viewd_supers::HostMirror initialise_supers_on_host(const GenSuperdrop &gen,
                                                   const viewd_supers supers);

/**
 * @brief Check if superdroplets initialisation is complete and sorted.
 *
 * This function checks if the initialisation of supers view is complete and if the super-droplets
 * are sorted by ascending superdroplet Gridbox index. If the initialisation is incomplete
 * or the superdroplets are not sorted, it throws an exception with an appropriate error message.
 *
 * @param supers The view of super-droplets in device memory.
 * @param size The expected number of super-droplets.
 *
 * @throws std::invalid_argument If the initialisation is incomplete or super-droplets
 * are not ordered correctly.
 */
void is_sdsinit_complete(const viewd_constsupers supers, const size_t size);

/**
 * @brief Print statement about initialised super-droplets.
 *
 * This function prints information about each superdroplet, including its ID, Gridbox index,
 * spatial coordinates, and attributes.
 *
 * @param supers The view of super-droplets in device memory.
 */
void print_supers(const viewd_constsupers supers);

/**
 * @brief Create a view of super-droplets in (device) memory.
 *
 * This function creates an ordered view of superdrops in device memory, where the number
 * of superdrops is specified by the parameter `totnsupers`. The superdrops are
 * ordered by the gridbox indexes and generated using a generator which uses
 * the initial conditions referenced by the `SuperdropInitConds` type.
 *
 * @tparam SuperdropInitConds The type of the super-droplets' initial conditions data.
 * @param sdic The instance of the super-droplets' initial conditions data.
 * @return A view of super-droplets in device memory.
 */
template <typename SuperdropInitConds>
viewd_supers create_supers(const SuperdropInitConds &sdic) {
  // Log message and create superdrops using the initial conditions
  std::cout << "\n--- create superdrops ---\ninitialising\n";
  viewd_supers supers(initialise_supers(sdic));

  // Log message and sort the view of superdrops
  std::cout << "sorting\n";
  supers = sort_supers(supers);

  // Log message and perform checks on the initialisation of superdrops
  std::cout << "checking initialisation\n";
  is_sdsinit_complete(supers, sdic.fetch_data_size());

  // Print information about the created superdrops
  print_supers(supers);

  // Log message indicating the successful creation of superdrops
  std::cout << "--- create superdrops: success ---\n";

  return supers;
}

/**
 * @brief Return an initialised view of superdrops in device memory.
 *
 * This function initialises a view of superdrops in device memory by creating
 * a view on the device and copying a host mirror view that is initialised using
 * the `SuperdropInitConds` instance.
 *
 * @tparam SuperdropInitConds The type of the super-droplets' initial conditions data.
 * @param sdic The instance of the super-droplets' initial conditions data.
 * @return A view of superdrops in device memory.
 */
template <typename SuperdropInitConds>
viewd_supers initialise_supers(const SuperdropInitConds &sdic) {
  // create superdrops view on device
  viewd_supers supers("supers", sdic.get_totnsupers());

  // initialise a mirror of superdrops view on host
  const GenSuperdrop gen(sdic);
  auto h_supers = initialise_supers_on_host(gen, supers);

  // Copy host view to device (h_supers to supers)
  Kokkos::deep_copy(supers, h_supers);

  return supers;
}

/**
 * @brief Return a mirror view of superdrops on host memory.
 *
 * This function initialises a mirror view of superdrops on host memory, using the super-droplet
 * generator instance `GenSuperdrop` to generate the kk'th super-droplet with their initial Gridbox
 * index, spatial coordinates, and attributes. Kokkos::parallel_for is used to perform parallel
 * initialisation of the mirror view, where each superdrop is generated by the provided `gen`
 * function object.
 *
 * The equivalent serial version of the Kokkos::parallel_for([...]) is:
 * @code
 * for (size_t kk(0); kk < totnsupers; ++kk)
 * {
 *  h_supers(kk) = gen(kk);
 * }
 * @endcode
 *
 * @param gen The generator function object used to initialize the superdrops.
 * @param supers The view of superdrops on device memory.
 * @return A mirror view of superdrops on host memory.
 */
viewd_supers::HostMirror initialise_supers_on_host(const GenSuperdrop &gen,
                                                   const viewd_supers supers) {
  const size_t totnsupers(supers.extent(0));

  // Create a mirror view of supers in case the original view is on device memory
  auto h_supers = Kokkos::create_mirror_view(supers);

  // Parallel initialisation of the mirror view
  Kokkos::parallel_for("initialise_supers_on_host", Kokkos::RangePolicy<HostSpace>(0, totnsupers),
                       [=](const size_t kk) { h_supers(kk) = gen(kk); });

  return h_supers;
}

#endif  // LIBS_RUNCLEO_CREATESUPERS_HPP_
