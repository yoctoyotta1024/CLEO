/*
 * ----- CLEO -----
 * File: superdrop.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 14th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for definition of a superdropet.
 * Equations referenced as (eqn [X.YY]) are from
 * "An Introduction To Clouds From The Microscale
 * to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 */

#ifndef SUPERDROP_HPP
#define SUPERDROP_HPP

#include <cmath>
#include <memory>
#include <utility>
#include <algorithm>
#include <string>
#include <stdexcept>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "./superdrop_ids.hpp"

namespace dlc = dimless_constants;

struct SoluteProperties
{
  const double rho_sol; // (dimensionless) density of solute in droplets
  const double mr_sol;   // (dimensionless) Mr of solute
  const double ionic;   // degree ionic dissociation (van't Hoff factor)

  /* A Kokkos requirement for use of struct in (dual)View (such as a
  Kokkos::vector) is that default constructor and destructor exist */
  KOKKOS_INLINE_FUNCTION
  SoluteProperties()
      : rho_sol(dlc::Rho_sol),
        mr_sol(dlc::Mr_sol),
        ionic(dlc::IONIC) {}

  KOKKOS_INLINE_FUNCTION ~SoluteProperties() = default;
};

struct SuperdropAttrs
/* attributes of a superdroplet*/
{
  double radius;                                  // radius of superdroplet
  double msol;                                    // mass of solute dissovled
  unsigned long long xi;                          // multiplicity of superdroplet
  std::shared_ptr<const SoluteProperties> solute; // reference to solute properties
};

class Superdrop
{
private:
  unsigned int sd_gbxindex; // matches index of gridbox the superdrop occupies
  double coord3;            // a 3rd spatial coordinate of superdroplet (z)
  double coord1;            // a 1st spatial coordinate of superdroplet (x)
  double coord2;            // a 2nd spatial coordinate of superdroplet (y)
  SuperdropAttrs attrs;     // attributes of a superdroplet

public:
  using IDType = IntID; // type of ID (from superdrop_ids.hpp) to identify superdrop via integer
  // using IDType = EmptyID; // empty type of ID (for non-existent superdrop identity)
  [[no_unique_address]] IDType id; // superdroplet (unique) identity

  KOKKOS_INLINE_FUNCTION Superdrop() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Superdrop() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  Superdrop(const unsigned int isd_gbxindex,
            const double icoord3,
            const double icoord1,
            const double icoord2,
            const double iradius,
            const double imsol,
            const unsigned long long ixi,
            const std::shared_ptr<const SoluteProperties> isolute,
            const IDType iid)
      : sd_gbxindex(isd_gbxindex),
        coord3(icoord3), coord1(icoord1), coord2(icoord2),
        attrs({iradius, imsol, ixi, isolute}),
        id(iid) {}

  KOKKOS_INLINE_FUNCTION
  Superdrop(const unsigned int isd_gbxindex,
            const double icoord3,
            const double icoord1,
            const double icoord2,
            const SuperdropAttrs iattrs,
            const IDType iid)
      : sd_gbxindex(isd_gbxindex),
        coord3(icoord3), coord1(icoord1), coord2(icoord2),
        attrs(iattrs), id(iid) {} 

  KOKKOS_INLINE_FUNCTION unsigned int get_sd_gbxindex() { return sd_gbxindex; }
  KOKKOS_INLINE_FUNCTION double get_coord3() { return coord3; }
  KOKKOS_INLINE_FUNCTION double get_coord1() { return coord1; }
  KOKKOS_INLINE_FUNCTION double get_coord2() { return coord2; }
  KOKKOS_INLINE_FUNCTION double get_radius() { return attrs.radius; }
  KOKKOS_INLINE_FUNCTION double get_msol() { return attrs.msol; }
  KOKKOS_INLINE_FUNCTION double get_xi() { return attrs.xi; }
};

#endif // SUPERDROP_HPP