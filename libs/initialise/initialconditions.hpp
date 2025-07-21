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
#include <cstdint>
#include <utility>
#include <vector>

#include "../kokkosaliases.hpp"
#include "superdrops/superdrop.hpp"

/* appends vector b onto end of vector a*/
template <typename T>
inline std::vector<T> append_vector(const std::vector<T> a, const std::vector<T> b) {
  auto ab = a;
  ab.insert(ab.end(), b.begin(), b.end());
  return ab;
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
  std::vector<Superdrop::IDType> sdIds;

  InitSupersData operator+(const InitSupersData& other) const {
    auto solutes_ = solutes;
    auto sdgbxindexes_ = append_vector(sdgbxindexes, other.sdgbxindexes);
    auto coord3s_ = append_vector(coord3s, other.coord3s);
    auto coord1s_ = append_vector(coord1s, other.coord1s);
    auto coord2s_ = append_vector(coord2s, other.coord2s);
    auto radii_ = append_vector(radii, other.radii);
    auto msols_ = append_vector(msols, other.msols);
    auto xis_ = append_vector(xis, other.xis);
    auto sdIds_ = append_vector(sdIds, other.sdIds);

    return InitSupersData{solutes_, sdgbxindexes_, coord3s_, coord1s_, coord2s_,
                          radii_,   msols_,        xis_,     sdIds_};
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
concept InitialConditions = requires(IC ic) {
  { ic.initsupers.get_maxnsupers() } -> std::convertible_to<size_t>;
  { ic.initsupers.get_nspacedims() } -> std::convertible_to<unsigned int>;
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
