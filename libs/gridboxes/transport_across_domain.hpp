/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: transport_across_domain.hpp
 * Project: gridboxes
 * Created Date: Monday 3rd March 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept for types used by MoveSupersInDomain to transport superdroplets across the domain
 * (see movesupersindomain.hpp)
 */

#ifndef LIBS_GRIDBOXES_TRANSPORT_ACROSS_DOMAIN_HPP_
#define LIBS_GRIDBOXES_TRANSPORT_ACROSS_DOMAIN_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../kokkosaliases.hpp"
#include "gridboxes/supersindomain.hpp"

/*
 * concept for TransportAcrossDomain is all types that have correct signature
 * for the following functions (e.g. operator(...) )
 */
template <typename T, typename GbxMaps>
concept TransportAcrossDomain =
    requires(T t, const GbxMaps &gbxmaps, const viewd_gbx d_gbxs, SupersInDomain &allsupers) {
      { t(gbxmaps, d_gbxs, allsupers) } -> std::convertible_to<SupersInDomain>;
    };

#endif  // LIBS_GRIDBOXES_TRANSPORT_ACROSS_DOMAIN_HPP_
