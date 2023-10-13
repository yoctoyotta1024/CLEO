/*
 * ----- CLEO -----
 * File: superdrop.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
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
  const double rho_l;   // (dimensionless) density of liquid in droplets
  const double rho_sol; // (dimensionless) density of solute in droplets
  const double mrsol;   // (dimensionless) Mr of solute
  const double ionic;   // degree ionic dissociation (van't Hoff factor)

  /* A Kokkos requirement for use of struct in (dual)View (such as a
  Kokkos::vector) is that default constructor and destructor exist */
  KOKKOS_INLINE_FUNCTION
  SoluteProperties() : rho_l(dlc::Rho_l), rho_sol(dlc::Rho_sol),
                       mrsol(dlc::Mr_sol), ionic(dlc::IONIC) {}

  KOKKOS_INLINE_FUNCTION ~SoluteProperties() = default;
};

class Superdrop
{
private:
  std::shared_ptr<const SoluteProperties> solute; // reference to solute properties

public:
  using IDType = IntID; // type of ID (from superdrop_ids.hpp) to identify superdrop via integer
  // using IDType = EmptyID; // empty type of ID (for non-existent superdrop identity)

  unsigned int sd_gbxindex;        // matches index of gridbox the superdrop occupies
  unsigned long long eps;          // multiplicity of superdroplet
  double radius;                   // radius of superdroplet
  double m_sol;                    // mass of solute dissovled
  double coord3;                   // a 3rd spatial coordinate of superdroplet (z)
  double coord1;                   // a 1st spatial coordinate of superdroplet (x)
  double coord2;                   // a 2nd spatial coordinate of superdroplet (y)
  [[no_unique_address]] IDType id; // superdroplet (unique) identity

  KOKKOS_INLINE_FUNCTION
  Superdrop(const std::shared_ptr<const SoluteProperties> isolute,
            const unsigned int isd_gbxindex,
            const unsigned long long ieps, const double iradius,
            const double im_sol, const double icoord3,
            const double icoord1, const double icoord2,
            const IDType iid)
      : solute(isolute), sd_gbxindex(isd_gbxindex), eps(ieps), radius(iradius),
        m_sol(im_sol), coord3(icoord3), coord1(icoord1),
        coord2(icoord2), id(iid) {}

  KOKKOS_INLINE_FUNCTION Superdrop() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Superdrop() = default; // Kokkos requirement for a (dual)View
};

struct Superdrops
{
  // replace with kokkos vector 
};

#endif // SUPERDROP_HPP