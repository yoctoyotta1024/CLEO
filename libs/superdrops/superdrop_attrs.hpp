/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: superdrop_attrs.hpp
 * Project: superdrops
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 18th January 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for structs and functions for
 * attributes superdroplets (note this excludes
 * gridbox index, coordinates and unique ids)
 * e.g. for solute, radius, multiplicity etc.
 */

#ifndef LIBS_SUPERDROPS_SUPERDROP_ATTRS_HPP_
#define LIBS_SUPERDROPS_SUPERDROP_ATTRS_HPP_

#include <cassert>

#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalConstants.hpp>  // for pi

#include "../cleoconstants.hpp"

namespace dlc = dimless_constants;

/* pointer-like object for solute properties of superdrop */
struct SoluteProperties {
  /* (dimensionless) density of solute in droplets */
  KOKKOS_INLINE_FUNCTION constexpr double rho_sol() const { return dlc::Rho_sol; }

  /* (dimensionless) molecular mass of solute */
  KOKKOS_INLINE_FUNCTION constexpr double mr_sol() const { return dlc::Mr_sol; }

  /* degree ionic dissociation (van't Hoff factor) */
  KOKKOS_INLINE_FUNCTION constexpr double ionic() const { return dlc::IONIC; }
};

/* attributes of a superdroplet*/
struct SuperdropAttrs {
  SoluteProperties solute;  // pointer-like reference to properties of solute
  uint64_t xi;              // multiplicity of superdroplet
  double radius;            // radius of superdroplet
  double msol;              // mass of solute dissovled

  KOKKOS_INLINE_FUNCTION SuperdropAttrs() = default;   // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~SuperdropAttrs() = default;  // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  SuperdropAttrs(const SoluteProperties solute, const uint64_t xi, const double radius,
                 const double msol)
      : solute(solute), xi(xi), radius(radius), msol(msol) {}

  KOKKOS_INLINE_FUNCTION bool is_solute() const { return true; }  // true if solute is "allocated"
  KOKKOS_INLINE_FUNCTION auto get_solute() const { return solute; }
  KOKKOS_INLINE_FUNCTION auto get_rho_sol() const { return solute.rho_sol(); }
  KOKKOS_INLINE_FUNCTION auto get_mr_sol() const { return solute.mr_sol(); }
  KOKKOS_INLINE_FUNCTION auto get_ionic() const { return solute.ionic(); }

  KOKKOS_INLINE_FUNCTION
  void set_xi(const uint64_t i_xi) {
    assert((i_xi > 0) && "xi should not be less than 1");

    xi = i_xi;
  }

  /* see also change_radius which prevents
  drop radius < dry radius */
  KOKKOS_INLINE_FUNCTION
  void set_radius(const double i_radius) {
    assert((i_radius - dryradius() > -1e-12 / dlc::R0) &&
           "radius cannot be less than dry radius " && "(within 1e-6 micron tolerance)");

    radius = i_radius;
  }

  KOKKOS_INLINE_FUNCTION
  void set_msol(const double i_msol) { msol = i_msol; }

  /* returns total droplet mass = water + dry areosol  */
  KOKKOS_INLINE_FUNCTION double mass() const;

  /* calculate radius as if dry droplet, ie.
  radius if drop was entirely made of solute */
  KOKKOS_INLINE_FUNCTION double dryradius() const {
    constexpr double vconst = 3.0 / (4.0 * Kokkos::numbers::pi);
    const auto dryrcubed = double{vconst * msol / solute.rho_sol()};
    return Kokkos::pow(dryrcubed, 1.0 / 3.0);
  }

  KOKKOS_INLINE_FUNCTION double rcubed() const { return radius * radius * radius; }

  /* spherical volume of droplet calculated from its radius */
  KOKKOS_INLINE_FUNCTION double vol() const { return 4.0 / 3.0 * Kokkos::numbers::pi * rcubed(); }

  /* Update droplet radius to newr or dry_radius() and
  return resultant change in radius (delta_radius = newradius-radius).
  Prevents drops shrinking further once they are size of dry_radius(). */
  KOKKOS_INLINE_FUNCTION double change_radius(const double newr);
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

/* Update droplet radius to newr or dry_radius() and
return resultant change in radius (delta_radius = newradius-radius).
Prevents drops shrinking further once they are size of dry_radius(). */
KOKKOS_INLINE_FUNCTION
double SuperdropAttrs::change_radius(const double newr) {
  const auto oldradius = radius;

  /*  if droplets are dry, do not shrink further */
  const auto dryr = dryradius();
  radius = Kokkos::fmax(
      newr, dryr);  // Kokkos compatible equivalent to std::max() for floating point numbers

  /* return change in radius due to growth/shrinking of droplet */
  return radius - oldradius;
}

/* returns (dimensionless) total droplet mass = water + dry areosol  */
KOKKOS_INLINE_FUNCTION
double SuperdropAttrs::mass() const {
  constexpr double massconst(4.0 / 3.0 * Kokkos::numbers::pi * dlc::Rho_l);  // 4/3 * pi * density
  const auto density_factor = double{1.0 - dlc::Rho_l / solute.rho_sol()};   // to account for msol

  auto mass = double{msol * density_factor};  // mass contribution of solute
  mass += massconst * rcubed();

  return mass;
}

#endif  // LIBS_SUPERDROPS_SUPERDROP_ATTRS_HPP_
