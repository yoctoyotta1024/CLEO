/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: superdrop_attrs.hpp
 * Project: superdrops
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for structs for attributes of super-droplets (note this excludes a super-droplet's
 * gridbox index, coordinates and unique ID, and includes e.g. a super-droplet's solute, radius,
 * multiplicity etc.).
 */

#ifndef LIBS_SUPERDROPS_SUPERDROP_ATTRS_HPP_
#define LIBS_SUPERDROPS_SUPERDROP_ATTRS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalConstants.hpp>  // for pi
#include <cassert>
#include <cstdint>

#include "../cleoconstants.hpp"

namespace dlc = dimless_constants;

/**
 * @brief Struct representing the properties of solute in a super-droplet.
 */
struct SoluteProperties {
  /**
   * @brief Get the density of solute (dimensionless).
   *
   * @return Solute density.
   */
  KOKKOS_FUNCTION constexpr double rho_sol() const { return dlc::Rho_sol; }

  /**
   * @brief Get the molecular mass of solute (dimensionless).
   *
   * @return The molecular mass of solute.
   */
  KOKKOS_FUNCTION constexpr double mr_sol() const { return dlc::Mr_sol; }

  /**
   * @brief Get the degree ionic dissociation (van't Hoff factor) (dimensionless).
   *
   * @return The degree ionic dissociation.
   */
  KOKKOS_FUNCTION constexpr double ionic() const { return dlc::IONIC; }
};

/**
 * @brief Struct representing the attributes of a super-droplet.
 */
struct SuperdropAttrs {
  SoluteProperties solute;
  /**< instance of SoluteProperties for pointer-like reference to superdrop's solute properties. */
  uint64_t xi;   /**< Multiplicity of superdrop. */
  double radius; /**< Radius of superdrop. */
  double msol;   /**< Mass of solute dissolved in superdrop. */

  /**
   * @brief Default constructor requirement for use of SuperdropAttrs in Kokkos View
   */
  SuperdropAttrs() = default;

  /**
   * @brief Default destructor requirement for use of SuperdropAttrs in Kokkos View
   */
  ~SuperdropAttrs() = default;

  /**
   * @brief Constructor with parameters.
   *
   * @param solute The solute properties.
   * @param xi The multiplicity of superdroplet.
   * @param radius The radius of superdroplet.
   * @param msol The mass of solute dissolved.
   * @param allow_nans Boolean to allow non-real superdroplets
   */
  KOKKOS_FUNCTION
  SuperdropAttrs(const SoluteProperties solute, const uint64_t xi, const double radius,
                 const double msol, const bool allow_nans = false)
      : solute(solute), xi(xi), radius(radius), msol(msol) {
    if (!allow_nans) {
      // if un-real superdroplets are not allowed check setter functions succeed
      set_xi(xi);
      set_radius(radius);
      set_msol(msol);
    }
  }

  /**
   * @brief Check if solute is present.
   *
   * @return True
   */
  KOKKOS_FUNCTION bool is_solute() const { return true; }

  /**
   * @brief Get the solute properties.
   *
   * @return The solute properties.
   */
  KOKKOS_FUNCTION auto get_solute() const { return solute; }

  /**
   * @brief Get the density of solute.
   *
   * @return The solute density.
   */
  KOKKOS_FUNCTION auto get_rho_sol() const { return solute.rho_sol(); }

  /**
   * @brief Get the molecular mass of the solute.
   *
   * @return The molecular mass of the solute.
   */
  KOKKOS_FUNCTION auto get_mr_sol() const { return solute.mr_sol(); }

  /**
   * @brief Get the degree ionic dissociation (van't Hoff factor).
   *
   * @return The ionic dissociation (van't Hoff factor).
   */
  KOKKOS_FUNCTION auto get_ionic() const { return solute.ionic(); }

  /**
   * @brief Set the multiplicity of superdroplet.
   *
   * Set the multiplicity 'xi' of the superdroplet with assert that new xi >= 1.
   *
   * @param i_xi The multiplicity to set.
   */
  KOKKOS_FUNCTION
  void set_xi(const uint64_t i_xi) {
    assert((i_xi > 0) && "xi should not be less than 1");

    xi = i_xi;
  }

  /**
   * @brief Set the radius of the super-droplet.
   *
   * This function sets the value of the super-droplet's radius to the specified value
   * with assert that new radius >= dry radius within 10^(-6) micron tolerance.
   *
   * _Note:_ See also change_radius which limits super-droplet radius to its dry radius.
   *
   * @param i_radius The value to set for radius.
   */
  KOKKOS_FUNCTION
  void set_radius(const double i_radius) {
    assert((i_radius - dryradius() > -1e-12 / dlc::R0) &&
           "radius cannot be less than dry radius (within 1e-6 micron tolerance)");

    radius = i_radius;
  }

  /**
   * @brief Sets the value of the super-droplet's mass of solute.
   *
   * This function sets the value of the super-droplet's solute mass to the specified value.
   *
   * @param i_msol The value to set for msol.
   */
  KOKKOS_FUNCTION
  void set_msol(const double i_msol) { msol = i_msol; }

  /**
   * @brief Get the total droplet mass.
   *
   * Calculates and returns total droplet mass = water + dry areosol.
   *
   * @return The total droplet mass.
   */
  KOKKOS_FUNCTION double mass() const;

  /**
   * @brief Get the mass of the droplet excluding its solute.
   *
   * @return mass of the super-droplet - mass of solute
   */
  KOKKOS_INLINE_FUNCTION double condensate_mass() const {
    auto m_cond = mass() - msol;
    assert((m_cond > -0.0001 * msol) &&
           "condensate mass cannot be less than 0.0 (within 0.0001 of dry mass tolerance)");

    m_cond = Kokkos::fmax(0.0, m_cond);  // Kokkos version of std::max() for floats (gpu compatible)

    return m_cond;
  }

  /**
   * @brief Get the dry radius of droplet.
   *
   * Calculate radius as if droplet is dry, ie. radius if droplet was only made of its solute mass
   *
   * @return The dry radius of droplet.
   */
  KOKKOS_FUNCTION double dryradius() const {
    constexpr double vconst = 3.0 / (4.0 * Kokkos::numbers::pi);
    const auto dryrcubed = double{vconst * msol / solute.rho_sol()};
    return Kokkos::pow(dryrcubed, 1.0 / 3.0);
  }

  /**
   * @brief Get the cubed radius of droplet.
   *
   * @return The cubed radius of droplet.
   */
  KOKKOS_FUNCTION double rcubed() const { return radius * radius * radius; }

  /**
   * @brief Get the spherical volume of droplet.
   *
   * Volume as if droplet is sphere given its radius.
   *
   * @return The spherical volume of droplet.
   */
  KOKKOS_FUNCTION double vol() const { return 4.0 / 3.0 * Kokkos::numbers::pi * rcubed(); }

  /**
   * @brief Change the radius of droplet.
   *
   * Update droplet radius to larger out of new radius 'newr' or dry radius and return the
   * resultant change in radius = new radius - old radius. Prevents drops shrinking further once
   * they are size of dry radius.
   *
   * @param newr The new radius to set.
   * @return The change in radius.
   */
  KOKKOS_FUNCTION double change_radius(const double newr);
};

#endif  // LIBS_SUPERDROPS_SUPERDROP_ATTRS_HPP_
