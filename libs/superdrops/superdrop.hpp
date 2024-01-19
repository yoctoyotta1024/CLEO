/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: superdrop.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 30th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for definition of a superdropet.
 * Equations referenced as (eqn [X.YY]) are from
 * "An Introduction To Clouds From The Microscale
 * to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 */

#ifndef SUPERDROP_HPP
#define SUPERDROP_HPP

#include "./superdrop_ids.hpp"
#include "./superdrop_attrs.hpp"

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
  [[no_unique_address]] IDType sdId; // superdroplet (unique) identity

  KOKKOS_INLINE_FUNCTION Superdrop() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Superdrop() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  Superdrop(const unsigned int i_sdgbxindex,
            const double i_coord3,
            const double i_coord1,
            const double i_coord2,
            const SuperdropAttrs i_attrs,
            const IDType i_sdId)
      : sdgbxindex(i_sdgbxindex),
        coord3(i_coord3),
        coord1(i_coord1),
        coord2(i_coord2),
        attrs(i_attrs),
        sdId(i_sdId) {}

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

  KOKKOS_INLINE_FUNCTION double mass() const { return attrs.mass(); }
  KOKKOS_INLINE_FUNCTION double vol() const { return attrs.vol(); }

  KOKKOS_INLINE_FUNCTION double rcubed() const { return attrs.rcubed(); }

  KOKKOS_INLINE_FUNCTION
  void set_xi(const unsigned long long i_xi)
  {
    attrs.set_xi(i_xi);
  }

  KOKKOS_INLINE_FUNCTION
  void set_radius(const double i_radius)
  /* see also change_radius which prevents
  drop radius < dry radius */
  {
    attrs.set_radius(i_radius);
  }

  KOKKOS_INLINE_FUNCTION
  void set_msol(const double i_msol)
  {
    attrs.set_msol(i_msol);
  }

  KOKKOS_INLINE_FUNCTION
  double change_radius(const double newr)
  {
    return attrs.change_radius(newr);
  }

  KOKKOS_INLINE_FUNCTION
  void set_sdgbxindex(const unsigned int i_sdgbxindex)
  {
    sdgbxindex = i_sdgbxindex;
  }

  KOKKOS_INLINE_FUNCTION
  void set_coord3(const double i_coord3)
  {
    coord3 = i_coord3;
  }

  KOKKOS_INLINE_FUNCTION
  void set_coord1(const double i_coord1)
  {
    coord1 = i_coord1;
  }

  KOKKOS_INLINE_FUNCTION
  void set_coord2(const double i_coord2)
  {
    coord2 = i_coord2;
  }

  KOKKOS_INLINE_FUNCTION
  void increment_coords(const double delta3,
                        const double delta1,
                        const double delta2)
  {
    coord3 += delta3;
    coord1 += delta1;
    coord2 += delta2;
  }
};

#endif // SUPERDROP_HPP
