/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: boundary_conditions.hpp
 * Project: gridboxes
 * Created Date: Monday 3rd March 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept for types used by MoveSupersInDomain to apply boundary conditions
 * (see movesupersindomain.hpp)
 */

#ifndef LIBS_GRIDBOXES_BOUNDARY_CONDITIONS_HPP_
#define LIBS_GRIDBOXES_BOUNDARY_CONDITIONS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/supersindomain.hpp"

/*
 * concept for BoundaryConditions is all types that have correct signature
 * for the following functions (e.g. operator(...) )
 */
template <typename BCs, typename GbxMaps>
concept BoundaryConditions =
    requires(BCs b, const GbxMaps &gbxmaps, const viewd_gbx d_gbxs, SupersInDomain &allsupers) {
      { b.apply(gbxmaps, d_gbxs, allsupers) } -> std::convertible_to<SupersInDomain>;
    };

/*
 * struct satisfying BoundaryConditions concept for applying domain boundary conditions which
 * does nothing.
 */
struct NullBoundaryConditions {
  template <GridboxMaps GbxMaps>
  SupersInDomain apply(const GbxMaps &gbxmaps, viewd_gbx d_gbxs,
                       const SupersInDomain &allsupers) const {
    return allsupers;
  }
};

#endif  // LIBS_GRIDBOXES_BOUNDARY_CONDITIONS_HPP_
