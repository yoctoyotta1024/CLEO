/*
 * ----- CLEO -----
 * File: superdrop_attrs.hpp
 * Project: superdrops
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 26th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for structs and functions for
 * attributes superdroplets (note this excludes
 * gridbox index, coordinates and unique ids)
 * e.g. for solute, radius, multiplicity etc.
 */

#ifndef SUPERDROP_ATTRS_HPP
#define SUPERDROP_ATTRS_HPP

#include <cmath>
#include <math.h>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"

namespace dlc = dimless_constants;

struct SoluteProperties
/* pointer-like object for solute properties of superdrop */
{
  /* (dimensionless) density of solute in droplets */
  KOKKOS_INLINE_FUNCTION constexpr double rho_sol() const { return dlc::Rho_sol; }

  /* (dimensionless) molecular mass of solute */
  KOKKOS_INLINE_FUNCTION constexpr double mr_sol() const { return dlc::Mr_sol; }

  /* degree ionic dissociation (van't Hoff factor) */
  KOKKOS_INLINE_FUNCTION constexpr double ionic() const { return dlc::IONIC; }
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

  KOKKOS_INLINE_FUNCTION bool is_solute() const { return true; } // true if solute is "allocated"
  KOKKOS_INLINE_FUNCTION auto get_solute() const { return solute; }
  KOKKOS_INLINE_FUNCTION auto get_rho_sol() const { return solute.rho_sol(); }
  KOKKOS_INLINE_FUNCTION auto get_mr_sol() const { return solute.mr_sol(); }
  KOKKOS_INLINE_FUNCTION auto get_ionic() const { return solute.ionic(); }

  KOKKOS_FUNCTION double mass() const;
  /* returns total droplet mass = water + dry areosol  */

  KOKKOS_INLINE_FUNCTION double dryradius() const
  /* calculate radius as if dry droplet, ie.
  radius if drop was entirely made of solute */
  {
    constexpr double vconst = 3.0 / (4.0 * M_PI);
    const double dryrcubed = vconst * msol / solute.rho_sol();
    return std::pow(dryrcubed, 1.0 / 3.0);
  }

  KOKKOS_FUNCTION double change_radius(const double newr);
  /* Update droplet radius to newr or dry_radius() and
  return resultant change in radius (delta_radius = newradius-radius).
  Prevents drops shrinking further once they are size of dry_radius(). */
};

#endif // SUPERDROP_ATTRS_HPP