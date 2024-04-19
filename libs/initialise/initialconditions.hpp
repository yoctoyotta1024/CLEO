/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: initialconditions.hpp
 * Project: initialise
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 19th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * structures and concept for generator of initial conditions for super-droplets and Gridboxes
 * in RunCLEO.
 */

#ifndef LIBS_INITIALISE_INITIALCONDITIONS_HPP_
#define LIBS_INITIALISE_INITIALCONDITIONS_HPP_

#include <array>
#include <concepts>
#include <utility>
#include <vector>

#include "../kokkosaliases.hpp"
#include "superdrops/superdrop_attrs.hpp"

template <typename T>
inline std::vector<T> join_vectors(const std::vector<T> a, const std::vector<T> b) {
  auto ab = a;
  ab.insert(ab.end(), b.begin(), b.end());
  return ab
}

/* struct required to generate initial super-droplets (see GenSuperdrop) */
struct InitSupersData {
  std::array<SoluteProperties, 1> solutes;
  std::vector<unsigned int> sdgbxindexes;
  std::vector<double> coord3s;
  std::vector<double> coord1s;
  std::vector<double> coord2s;
  std::vector<double> radii;
  std::vector<double> msols;
  std::vector<uint64_t> xis;

  InitSupersData operator+(const InitSupersData& other) const {
    auto solutes_ = solutes;
    auto sdgbxindexes_ = join_vectors(sdgbxindexes, other.sdgbxindexes);
    auto coord3s_ = join_vectors(coord3s, other.coord3s);
    auto coord1s_ = join_vectors(coord1s, other.coord1s);
    auto coord2s_ = join_vectors(coord2s, other.coord2s);
    auto radii_ = join_vectors(radii, other.radii);
    auto msols_ = join_vectors(msols, other.msols);
    auto xis_ = join_vectors(xis, other.xis);

    auto isd = InitSupersData{sdgbxindexes_, coord3s_, coord1s_, coord2s_, radii_, msols_, xis_};

    return isd;
  }
};

/**
 * @concept InitialConditions
 * Concept representing types that provide initial conditions used by operator
 * call of RunCLEO.
 *
 * A type satisfies the InitialConditions concept if it provides the following structures:
 * - `initsupers`: Struct that can call functions related to initialising super-droplets.
 * - `initgbxs`: Struct that can call functions related to initialising Gridboxes.
 *
 * @tparam IC The type to check against the InitialConditions concept.
 */
template <typename IC>
concept InitialConditions =
    requires(IC ic, unsigned int t, const viewh_constgbx h_gbxs, InitSupersData initdata) {
      { ic.initsupers.get_nspacedims() } -> std::convertible_to<unsigned int>;
      { ic.initsupers.get_maxnsupers() } -> std::convertible_to<size_t>;
      { ic.initsupers.fetch_data() } -> std::same_as<InitSupersData>;

      { ic.initgbxs.get_ngbxs() } -> std::convertible_to<size_t>;
      { ic.initgbxs.press() } -> std::convertible_to<std::vector<double>>;
      { ic.initgbxs.temp() } -> std::convertible_to<std::vector<double>>;
      { ic.initgbxs.qvap() } -> std::convertible_to<std::vector<double>>;
      { ic.initgbxs.qcond() } -> std::convertible_to<std::vector<double>>;
      { ic.initgbxs.wvel() } -> std::convertible_to<std::vector<std::pair<double, double>>>;
      { ic.initgbxs.uvel() } -> std::convertible_to<std::vector<std::pair<double, double>>>;
      { ic.initgbxs.vvel() } -> std::convertible_to<std::vector<std::pair<double, double>>>;
    };

/* helpful struct satisyfing InitialConditions concept for functions to generate
initial conditions for CLEO */
template <typename SuperdropInitConds, typename GbxInitConds>
struct InitConds {
  SuperdropInitConds initsupers;  // initial conditions for creating superdroplets
  GbxInitConds initgbxs;          // initial conditions for creating gridboxes

  InitConds(const SuperdropInitConds initsupers, const GbxInitConds initgbxs)
      : initsupers(initsupers), initgbxs(initgbxs) {}
};

#endif  // LIBS_INITIALISE_INITIALCONDITIONS_HPP_
