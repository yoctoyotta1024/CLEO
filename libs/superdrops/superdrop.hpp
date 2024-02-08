/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: superdrop.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 8th February 2024
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

#ifndef LIBS_SUPERDROPS_SUPERDROP_HPP_
#define LIBS_SUPERDROPS_SUPERDROP_HPP_

#include "./superdrop_attrs.hpp"
#include "./superdrop_ids.hpp"

/**
 * @brief Class representing a super-droplet (synonyms: superdroplet, superdrop, SD).
 *
 * This class defines properties and operations of a super-droplet
 * (synonyms: superdroplet, superdrop, SD).
 */
class Superdrop {
 private:
  unsigned int sdgbxindex;  /**< Index of the gridbox the superdrop occupies */
  double coord3;            /**< 3rd spatial coordinate of the superdroplet (z) */
  double coord1;            /**< 1st spatial coordinate of the superdroplet (x) */
  double coord2;            /**< 2nd spatial coordinate of the superdroplet (y) */
  SuperdropAttrs attrs;     /**< Attributes of the superdroplet */

 public:
  using IDType = IntID;     /**< Type of ID to identify superdrop via 8 byte integer */
  // using IDType = EmptyID; // Type of ID (for non-existent superdrop identity)
  [[no_unique_address]] IDType sdId;  // superdroplet (unique) identity

  /**
   * @brief Default constructor requirement for use of Superdrop in Kokkos View
   */
  KOKKOS_INLINE_FUNCTION Superdrop() = default;

  /**
   * @brief Default destructor requirement for use of Superdrop in Kokkos View
   */
  KOKKOS_INLINE_FUNCTION ~Superdrop() = default;

  /**
   * @brief Parameterized constructor.
   *
   * @param i_sdgbxindex Index of the gridbox the superdrop occupies.
   * @param i_coord3 3rd spatial coordinate of the superdroplet.
   * @param i_coord1 1st spatial coordinate of the superdroplet.
   * @param i_coord2 2nd spatial coordinate of the superdroplet.
   * @param i_attrs Attributes of the superdroplet.
   * @param i_sdId Unique ID of the superdroplet.
   */
  KOKKOS_INLINE_FUNCTION
  Superdrop(const unsigned int i_sdgbxindex, const double i_coord3, const double i_coord1,
            const double i_coord2, const SuperdropAttrs i_attrs, const IDType i_sdId)
      : sdgbxindex(i_sdgbxindex),
        coord3(i_coord3),
        coord1(i_coord1),
        coord2(i_coord2),
        attrs(i_attrs),
        sdId(i_sdId) {}

  /**
   * @brief Get the index of the Gridbox the superdrop currently occupies.
   *
   * @return Gridbox Index.
   */
  KOKKOS_INLINE_FUNCTION auto get_sdgbxindex() const { return sdgbxindex; }

  /**
   * @brief Get the 3rd spatial coordinate of the superdroplet.
   *
   * @return 3rd spatial coordinate.
   */
  KOKKOS_INLINE_FUNCTION auto get_coord3() const { return coord3; }

  /**
   * @brief Get the 1st spatial coordinate of the superdroplet.
   *
   * @return 1st spatial coordinate.
   */
  KOKKOS_INLINE_FUNCTION auto get_coord1() const { return coord1; }

  /**
   * @brief Get the 2nd spatial coordinate of the superdroplet.
   *
   * @return 2nd spatial coordinate.
   */
  KOKKOS_INLINE_FUNCTION auto get_coord2() const { return coord2; }

  /**
   * @brief Returns 'true' if super-droplet has solute.
   *
   * @return boolean for existance of solute.
   */
  KOKKOS_INLINE_FUNCTION auto is_solute() const { return attrs.is_solute(); }

  /**
   * @brief Get the solute of the superdroplet.
   *
   * @return super-droplet's solute
   */
  KOKKOS_INLINE_FUNCTION auto get_solute() const { return attrs.get_solute(); }

  /**
   * @brief Get the density of the super-droplet's solute.
   *
   * @return density of the super-droplet's solute.
   */
  KOKKOS_INLINE_FUNCTION auto get_rho_sol() const { return attrs.get_rho_sol(); }

  /**
   * @brief Get the molecular mass of the super-droplet's solute.
   *
   * @return molecular mass of the super-droplet's solute.
   */
  KOKKOS_INLINE_FUNCTION auto get_mr_sol() const { return attrs.get_mr_sol(); }

  /**
   * @brief Get the van't Hoff ionic factor of the super-droplet's solute.
   *
   * @return van't Hoff ionic factor of the super-droplet's solute.
   */
  KOKKOS_INLINE_FUNCTION auto get_ionic() const { return attrs.get_ionic(); }

  /**
   * @brief Get the multiplicity 'xi' of the super-droplet.
   *
   * @return multiplicity 'xi' of the super-droplet.
   */
  KOKKOS_INLINE_FUNCTION auto get_xi() const { return attrs.xi; }

  /**
   * @brief Get the spherical radius of the super-droplet.
   *
   * @return spherical radius of the super-droplet.
   */
  KOKKOS_INLINE_FUNCTION auto get_radius() const { return attrs.radius; }

  /**
   * @brief Get the mass of solute dissolved in the super-droplet.
   *
   * @return mass of solute in the super-droplet.
   */
  KOKKOS_INLINE_FUNCTION auto get_msol() const { return attrs.msol; }

  /**
   * @brief Get the mass of the super-droplet.
   *
   * @return mass of the super-droplet.
   */
  KOKKOS_INLINE_FUNCTION double mass() const { return attrs.mass(); }

  /**
   * @brief Get the volume of the super-droplet.
   *
   * @return volume of the super-droplet.
   */
  KOKKOS_INLINE_FUNCTION double vol() const { return attrs.vol(); }

  /**
   * @brief Get the radius of the super-droplet cubed.
   *
   * @return radius of the super-droplet cubed.
   */
  KOKKOS_INLINE_FUNCTION double rcubed() const { return attrs.rcubed(); }

  KOKKOS_INLINE_FUNCTION
  void set_xi(const uint64_t i_xi) { attrs.set_xi(i_xi); }

  /* see also change_radius which prevents
  drop radius < dry radius */
  KOKKOS_INLINE_FUNCTION
  void set_radius(const double i_radius) { attrs.set_radius(i_radius); }

  KOKKOS_INLINE_FUNCTION
  void set_msol(const double i_msol) { attrs.set_msol(i_msol); }

  KOKKOS_INLINE_FUNCTION
  double change_radius(const double newr) { return attrs.change_radius(newr); }

  KOKKOS_INLINE_FUNCTION
  void set_sdgbxindex(const unsigned int i_sdgbxindex) { sdgbxindex = i_sdgbxindex; }

  KOKKOS_INLINE_FUNCTION
  void set_coord3(const double i_coord3) { coord3 = i_coord3; }

  KOKKOS_INLINE_FUNCTION
  void set_coord1(const double i_coord1) { coord1 = i_coord1; }

  KOKKOS_INLINE_FUNCTION
  void set_coord2(const double i_coord2) { coord2 = i_coord2; }

  KOKKOS_INLINE_FUNCTION
  void increment_coords(const double delta3, const double delta1, const double delta2) {
    coord3 += delta3;
    coord1 += delta1;
    coord2 += delta2;
  }
};

#endif  // LIBS_SUPERDROPS_SUPERDROP_HPP_
