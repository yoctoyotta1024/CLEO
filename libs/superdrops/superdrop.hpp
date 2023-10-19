/*
 * ----- CLEO -----
 * File: superdrop.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 19th October 2023
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

#include <memory>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "./superdrop_ids.hpp"

namespace dlc = dimless_constants;

struct SoluteProperties
{
  double rho_sol; // (dimensionless) density of solute in droplets
  double mr_sol;  // (dimensionless) Mr of solute
  double ionic;   // degree ionic dissociation (van't Hoff factor)

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
  using SolutePtr = Kokkos::View<const SoluteProperties[1]>; // pointer-like type to an instance of Solute Properties
  SolutePtr soluteptr;
  unsigned long long xi;                            // multiplicity of superdroplet
  double radius;                                    // radius of superdroplet
  double msol;                                      // mass of solute dissovled

  KOKKOS_INLINE_FUNCTION SuperdropAttrs() = default; // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~SuperdropAttrs() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  SuperdropAttrs(const SolutePtr soluteptr,
                 const unsigned long long xi,
                 const double radius,
                 const double msol)
      : soluteptr(soluteptr),
        xi(xi),
        radius(radius),
        msol(msol) {}
};

class Superdrop
{
private:
  unsigned int sdgbxindex; // matches index of gridbox the superdrop occupies
  double coord3;           // a 3rd spatial coordinate of superdroplet (z)
  double coord1;           // a 1st spatial coordinate of superdroplet (x)
  double coord2;           // a 2nd spatial coordinate of superdroplet (y)
  SuperdropAttrs attrs;    // attributes of a superdroplet

public:
  using IDType = IntID; // type of ID (from superdrop_ids.hpp) to identify superdrop via integer
  // using IDType = EmptyID; // empty type of ID (for non-existent superdrop identity)
  [[no_unique_address]] IDType id; // superdroplet (unique) identity

  KOKKOS_INLINE_FUNCTION Superdrop() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Superdrop() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  Superdrop(const unsigned int isdgbxindex,
            const double icoord3,
            const double icoord1,
            const double icoord2,
            const SuperdropAttrs iattrs,
            const IDType iid)
      : sdgbxindex(isdgbxindex),
        coord3(icoord3),
        coord1(icoord1),
        coord2(icoord2),
        attrs(iattrs),
        id(iid) {}

  KOKKOS_INLINE_FUNCTION auto get_sdgbxindex() const { return sdgbxindex; }
  KOKKOS_INLINE_FUNCTION auto get_coord3() const { return coord3; }
  KOKKOS_INLINE_FUNCTION auto get_coord1() const { return coord1; }
  KOKKOS_INLINE_FUNCTION auto get_coord2() const { return coord2; }

  KOKKOS_INLINE_FUNCTION auto is_soluteptr() const { return attrs.soluteptr.is_allocated(); }
  KOKKOS_INLINE_FUNCTION auto get_solute() const { return attrs.soluteptr(0); }
  KOKKOS_INLINE_FUNCTION auto get_rho_sol() const { return attrs.soluteptr(0).rho_sol; }
  KOKKOS_INLINE_FUNCTION auto get_mr_sol() const { return attrs.soluteptr(0).mr_sol; }
  KOKKOS_INLINE_FUNCTION auto get_ionic() const { return attrs.soluteptr(0).ionic; }

  KOKKOS_INLINE_FUNCTION auto get_xi() const { return attrs.xi; }
  KOKKOS_INLINE_FUNCTION auto get_radius() const { return attrs.radius; }
  KOKKOS_INLINE_FUNCTION auto get_msol() const { return attrs.msol; }
};

#endif // SUPERDROP_HPP