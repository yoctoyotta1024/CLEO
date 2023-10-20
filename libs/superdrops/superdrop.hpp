/*
 * ----- CLEO -----
 * File: superdrop.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
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
  KOKKOS_INLINE_FUNCTION
  double rho_sol() const
  /* (dimensionless) density of solute in droplets */
  {
    return dlc::Rho_sol;
  }

  KOKKOS_INLINE_FUNCTION
  double mr_sol() const
  /* (dimensionless) Mr of solute */
  {
    return dlc::Mr_sol;
  }

  KOKKOS_INLINE_FUNCTION
  double ionic() const
  /* degree ionic dissociation (van't Hoff factor) */
  {
    return dlc::IONIC;
  }
};

struct SuperdropAttrs
/* attributes of a superdroplet*/
{
  SoluteProperties solute; // pointer-like reference to properties of solute
  unsigned long long xi;   // multiplicity of superdroplet
  double radius;           // radius of superdroplet
  double msol;             // mass of solute dissovled

  KOKKOS_INLINE_FUNCTION SuperdropAttrs() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~SuperdropAttrs() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  SuperdropAttrs(const SoluteProperties solute,
                 const unsigned long long xi,
                 const double radius,
                 const double msol)
      : solute(solute),
        xi(xi),
        radius(radius),
        msol(msol) {}

  KOKKOS_INLINE_FUNCTION bool is_solute() const { return true; }
  KOKKOS_INLINE_FUNCTION auto get_solute() const { return solute; }
  KOKKOS_INLINE_FUNCTION auto get_rho_sol() const { return solute.rho_sol(); }
  KOKKOS_INLINE_FUNCTION auto get_mr_sol() const { return solute.mr_sol(); }
  KOKKOS_INLINE_FUNCTION auto get_ionic() const { return solute.ionic(); }
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

  KOKKOS_INLINE_FUNCTION auto is_solute() const { return attrs.is_solute(); }
  KOKKOS_INLINE_FUNCTION auto get_solute() const { return attrs.get_solute(); }
  KOKKOS_INLINE_FUNCTION auto get_rho_sol() const { return attrs.get_rho_sol(); }
  KOKKOS_INLINE_FUNCTION auto get_mr_sol() const { return attrs.get_mr_sol(); }
  KOKKOS_INLINE_FUNCTION auto get_ionic() const { return attrs.get_ionic(); }

  KOKKOS_INLINE_FUNCTION auto get_xi() const { return attrs.xi; }
  KOKKOS_INLINE_FUNCTION auto get_radius() const { return attrs.radius; }
  KOKKOS_INLINE_FUNCTION auto get_msol() const { return attrs.msol; }
};

#endif // SUPERDROP_HPP