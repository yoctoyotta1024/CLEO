// Author: Clara Bayley
// File: superdrop.hpp
/* Header file for superdrops class and
  structure for a superdrop with an
  associated gridbox. Equations referenced
  as (eqn [X.YY]) are from "An Introduction
  To Clouds From The Microscale to Climate"
  by Lohmann, Luond and Mahrt, 1st edition. */

#ifndef SUPERDROP_HPP
#define SUPERDROP_HPP

#include <cmath>
#include <memory>
#include <utility>
#include <algorithm>
#include <string>
#include <stdexcept>

#include <Kokkos_Core.hpp>

#include "../claras_SDconstants.hpp"
#include "./impliciteuler.hpp"
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

  KOKKOS_INLINE_FUNCTION double dry_radius() const
  /* calculate radius as if dry droplet, ie.
  radius if drop was entirely made of solute */
  {
    return pow(3.0 * m_sol / (4.0 * M_PI * solute->rho_sol), 1.0 / 3.0);
  }

  KOKKOS_FUNCTION double rhoeff() const;
  /* calculates effective density of droplet
  so mass_droplet = m = 4/3*pi*r^3 * rhoeff */

  // KOKKOS_INLINE_FUNCTION double mass_liq() const;
  // /* calculate mass of only water in droplet */

public:
  using IDType = IntID; // type of ID (from superdrop_ids.hpp) to identify superdrop via integer
  //using IDType = EmptyID; // empty type of ID (for non-existent superdrop identity)

  unsigned long long eps;    // multiplicity of superdroplet
  double radius; // radius of superdroplet
  double m_sol;  // mass of solute dissovled
  double coord3; // a 3rd spatial coordinate of superdroplet (z)
  double coord1; // a 1st spatial coordinate of superdroplet (x)
  double coord2; // a 2nd spatial coordinate of superdroplet (y)
  [[no_unique_address]] IDType id; // superdroplet (unique) identity

  KOKKOS_FUNCTION Superdrop(const std::shared_ptr<const SoluteProperties> isolute,
                            const unsigned long long ieps, const double iradius,
                            const double im_sol, const double icoord3,
                            const double icoord1, const double icoord2,
                            const IDType iid)
      : solute(isolute), eps(ieps), radius(iradius),
        m_sol(im_sol), coord3(icoord3), coord1(icoord1),
        coord2(icoord2), id(iid) {}

  KOKKOS_INLINE_FUNCTION Superdrop() = default; // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Superdrop() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION double vol() const
  /* spherical volume of droplet calculated from its radius */
  {
    return 4.0 / 3.0 * M_PI * pow(radius, 3.0);
  }

  KOKKOS_INLINE_FUNCTION double vol_liq() const
  /* volume of droplet excluding solute */
  {
    const double dryvol(m_sol / solute->rho_sol);

    return vol() - dryvol;
  }

  KOKKOS_FUNCTION double mass() const;
  /* calculate total mass of droplet
  mass = (water + dry areosol)  */

  KOKKOS_FUNCTION double equilibrium_wetradius(const double s_ratio,
                                               const double temp) const;
  /* Performs Newton Raphson root finding algorithm using functions in 
  WetRadius root finder struct to solve equation for equilibrium (wet) 
  radius of superdroplet at given relative humidity. Equilibrium radius 
  defined by radius when ODE from "An Introduction To Clouds...."
  (see note at top of file) eqn [7.28] = 0. */

  KOKKOS_FUNCTION double akohler_factor(double temp) const;
  /* calculate value of a in raoult factor (exp^(a/r))
  to account for effect of dissolved solute
  on radial growth of droplet. Using equations from
  "An Introduction To Clouds...." (see note at top of file) */

  KOKKOS_FUNCTION double bkohler_factor() const;
  /* calculate value of b in kelvin factor (1-b/r^3)
  to account for curvature on radial growth
  of droplet. Using equations from "An Introduction
  To Clouds...." (see note at top of file) */

  KOKKOS_FUNCTION double change_radius(const double newradius);
  /* Update droplet radius to newradius or dry_radius() and
  return resultant change in radius. Prevents drops shrinking
  further once they are size of dry_radius(). */

  KOKKOS_INLINE_FUNCTION auto get_solute() const
  {
    return solute;
  }

  KOKKOS_INLINE_FUNCTION double get_dry_radius() const
  /* calculate radius as if dry droplet, ie.
  radius if drop was entirely made of solute */
  {
    return dry_radius();
  }
};

struct WetRadius
{
  private:
    struct IterReturn
    /* struct used for returning boolean and double
    from iterate_rootfinding function */
    {
      bool do_iter;
      double ziter;
    };

    IterReturn iterate_rootfinding(double ziter, const double s_ratio,
                                   const double akoh, const double bkoh) const;
    /* iterate wetradius rootfinding algorithm. Performs 1 iteration of
    newton raphson root finding algorithm for obtaining the equilibrium
    wet radius of the condensation ODE at a given 
    relative humidity (s_ratio). ODE from "An Introduction To Clouds...." 
    (see note at top of file) eqn [7.28] */

    double wetradius_polynomial(const double ziter, const double s_ratio,
                                const double akoh, const double bkoh) const;
    /* returns value of (cubic) polynomial evaluted at ziter. Root of this
    polynomial is the value of the equilibrium (wet) radius of a superdroplet 
    at a given relative humidity (ie. s_ratio) derived from ODE in 
    "An Introduction To Clouds...." (see note at top of file) eqn [7.28] */                                    
    
    bool isnotconverged(const double new_ode, const double ode) const;
    /* boolean where True means criteria for ending newton raphson iterations
    has not yet been met. Criteria is standard local error test:
    |iteration - previous iteration| < RTOL * |iteration| + ATOL */
  
  public: 
    unsigned int maxiters;

    double get_wetradius(const double radius0, const double s_ratio,
                         const double akoh, const double bkoh) const;
    /* Iterate Newton Raphson root finding algorithm to
    return wet radius of a superdroplet in equilibrium
    with supersaturation s_ratio */
};

struct SuperdropWithGbxindex
/* Structure containing a superdroplet (SD) and the index/unique
identifier of the gridbox (GBx) it occupies ie. the identity of
the GBx the SD is associated with */
{
  unsigned int sd_gbxindex; // index/unique identifier of gridbox the superdrop occupies
  Superdrop superdrop;

  KOKKOS_INLINE_FUNCTION
  SuperdropWithGbxindex(const unsigned int isd_gbxindex, Superdrop isuperdrop)
      : sd_gbxindex(isd_gbxindex), superdrop(isuperdrop) {}
  
  KOKKOS_INLINE_FUNCTION SuperdropWithGbxindex() = default; // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~SuperdropWithGbxindex() = default; // Kokkos requirement for a (dual)View

};

#endif // SUPERDROP_HPP