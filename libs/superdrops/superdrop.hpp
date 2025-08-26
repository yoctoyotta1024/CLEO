/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: superdrop.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
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

#include <cstdint>

#include "superdrop_attrs.hpp"
#include "superdrop_ids.hpp"

/**
 * @brief Class representing a super-droplet (synonyms: superdroplet, superdrop, SD).
 *
 * This class defines properties and operations of a super-droplet
 * (synonyms: superdroplet, superdrop, SD).
 */
class Superdrop {
 private:
  unsigned int sdgbxindex; /**< Index of the gridbox the superdrop occupies */
  double coord3;           /**< 3rd spatial coordinate of the superdrop (vertical) */
  double coord1;           /**< 1st spatial coordinate of the superdrop (eastwards) */
  double coord2;           /**< 2nd spatial coordinate of the superdrop (northwards) */
  SuperdropAttrs attrs;    /**< instance of SuperdropAttrs for attributes of the super-droplet */

 public:
  // TODO(ALL): define superdrop IDType using macros(?)
  using IDType = IntID; /**< Type of ID to identify superdrop via 8 byte integer */
  // using IDType = EmptyID; /**< Type of ID for non-existent superdrop identity */
  [[no_unique_address]] IDType sdId;
  /**< instance of super-droplet identity of Superdrop::IDType */

  /**
   * @brief Default constructor requirement for use of Superdrop in Kokkos View
   */
  Superdrop() = default;

  /**
   * @brief Default destructor requirement for use of Superdrop in Kokkos View
   */
  ~Superdrop() = default;

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
   * @brief Returns 'true' if the super-droplet has solute.
   *
   * This function checks whether the super-droplet contains solute.
   *
   * @return true if the super-droplet has solute, false otherwise.
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
   * @brief Get the mass of the super-droplet excluding its solute.
   *
   * @return mass of the super-droplet - mass of solute
   */
  KOKKOS_INLINE_FUNCTION double condensate_mass() const { return attrs.condensate_mass(); }

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

  /**
   * @brief Set the multiplicity 'xi' of the super-droplet.
   *
   * This function sets the value of the super-droplet's multiplicity 'xi' to the specified value.
   *
   * @param i_xi The value to set for xi
   */
  KOKKOS_INLINE_FUNCTION
  void set_xi(const uint64_t i_xi) { attrs.set_xi(i_xi); }

  /**
   * @brief Set the radius of the super-droplet.
   *
   * This function sets the value of the super-droplet's radius to the specified value.
   *
   * _Note:_ See also change_radius which limits super-droplet radius to its dry radius.
   *
   * @param i_radius The value to set for radius.
   */
  KOKKOS_INLINE_FUNCTION
  void set_radius(const double i_radius) { attrs.set_radius(i_radius); }

  /**
   * @brief Sets the value of the super-droplet's mass of solute.
   *
   * This function sets the value of the super-droplet's solute mass to the specified value.
   *
   * @param i_msol The value to set for msol.
   */
  KOKKOS_INLINE_FUNCTION
  void set_msol(const double i_msol) { attrs.set_msol(i_msol); }

  /**
   * @brief Set the radius of the super-droplet no less than its dry radius.
   *
   * This function sets the value of the super-droplet's radius to the specified value
   * if its new radius >= its dry radius. Return is difference in radius accoding to super-droplet's
   * attributes struct.
   *
   * _Note:_ See also set_radius which allows super-droplet radius less than its dry radius
   *
   * @param newr The value to set the radius >= dry radius.
   * @return change in radius of the super-droplet.
   */
  KOKKOS_INLINE_FUNCTION
  double change_radius(const double newr) { return attrs.change_radius(newr); }

  /**
   * @brief Sets the value of the super-droplet's Gridbox index.
   *
   * This function sets the value of the super-droplet's sdgbxindex to the specified value.
   *
   * @param i_sdgbxindex The value to set for sdgbxindex.
   */
  KOKKOS_INLINE_FUNCTION
  void set_sdgbxindex(const unsigned int i_sdgbxindex) { sdgbxindex = i_sdgbxindex; }

  /**
   * @brief Sets the value of the 3rd coordinate.
   *
   * This function sets the value of coord3 of the super-droplet to the specified value.
   *
   * @param i_coord3 The value to set for coord3.
   */
  KOKKOS_INLINE_FUNCTION
  void set_coord3(const double i_coord3) { coord3 = i_coord3; }

  /**
   * @brief Sets the value of the 1st coordinate.
   *
   * This function sets the value of coord1 of the super-droplet to the specified value.
   *
   * @param i_coord1 The value to set for coord1.
   */
  KOKKOS_INLINE_FUNCTION
  void set_coord1(const double i_coord1) { coord1 = i_coord1; }

  /**
   * @brief Sets the value of the 2nd coordinate.
   *
   * This function sets the value of coord2 of the super-droplet to the specified value.
   *
   * @param i_coord2 The value to set for coord2.
   */
  KOKKOS_INLINE_FUNCTION
  void set_coord2(const double i_coord2) { coord2 = i_coord2; }

  /**
   * @brief Sets the values of the 3rd, 1st and 2nd coordinates.
   *
   * This function sets the coordinates of the super-droplet along each dimension.
   *
   * @param i_coord3 The value to set for coord3.
   * @param i_coord1 The value to set for coord1.
   * @param i_coord2 The value to set for coord2.
   */
  KOKKOS_INLINE_FUNCTION
  void set_coords(const double i_coord3, const double i_coord1, const double i_coord2) {
    coord3 = i_coord3;
    coord1 = i_coord1;
    coord2 = i_coord2;
  }

  /**
   * @brief Increments the coordinates by the specified deltas.
   *
   * This function increments the coordinates of the super-droplet by the specified deltas along
   * each dimension.
   *
   * @param delta3 The delta for the third coordinate (coord3).
   * @param delta1 The delta for the first coordinate (coord1).
   * @param delta2 The delta for the second coordinate (coord2).
   */
  KOKKOS_INLINE_FUNCTION
  void increment_coords(const double delta3, const double delta1, const double delta2) {
    coord3 += delta3;
    coord1 += delta1;
    coord2 += delta2;
  }

  void serialize_double_components(std::vector<double>::iterator target) const {
    *target++ = coord3;
    *target++ = coord1;
    *target++ = coord2;
    *target++ = attrs.radius;
    *target = attrs.msol;
  }

  void serialize_uint_components(std::vector<unsigned int>::iterator target) {
    *target++ = sdgbxindex;
    *target = static_cast<unsigned int>(sdId.get_value());  // TODO(ALL): don't do if using EmptyID
  }

  void serialize_uint64_components(std::vector<uint64_t>::iterator target) { *target = attrs.xi; }

  void deserialize_components(std::vector<unsigned int>::iterator uint_source,
                              std::vector<uint64_t>::iterator uint64_source,
                              std::vector<double>::iterator double_source) {
    sdgbxindex = *uint_source++;
    sdId.value = static_cast<size_t>(*uint_source);  // TODO(ALL): don't do if using EmptyID

    attrs.xi = *uint64_source;

    coord3 = *double_source++;
    coord1 = *double_source++;
    coord2 = *double_source++;
    attrs.radius = *double_source++;
    attrs.msol = *double_source;
  }
};

#endif  // LIBS_SUPERDROPS_SUPERDROP_HPP_
