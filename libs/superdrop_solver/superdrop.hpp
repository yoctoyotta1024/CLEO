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
#include <vector> 

#include "../claras_SDconstants.hpp"
#include "./impliciteuler.hpp"
#include "./superdrop_ids.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

struct SoluteProperties
{
  const double rho_l;   // (dimensionless) density of liquid in droplets
  const double rho_sol; // (dimensionless) density of solute in droplets
  const double mrsol;   // (dimensionless) Mr of solute
  const double ionic;   // degree ionic dissociation (van't Hoff factor)

  SoluteProperties() : rho_l(dlc::Rho_l), rho_sol(dlc::Rho_sol),
                       mrsol(dlc::Mr_sol), ionic(dlc::IONIC) {}
};

class Superdrop
{
private:
  std::shared_ptr<const SoluteProperties> solute; // reference to solute properties

  inline double dry_radius() const
  /* calculate radius as if dry droplet, ie.
  radius if drop was entirely made of solute */
  {
    return pow(3.0 * m_sol / (4.0 * M_PI * solute->rho_sol), 1.0 / 3.0);
  }

  double rhoeff() const;
  /* calculates effective density of droplet
  so mass_droplet = m = 4/3*pi*r^3 * rhoeff */

  // double m_liq() const;
  // /* calculate mass of only water in droplet */

public:
  using IDType = IntID; // type of ID (from superdrop_ids.hpp) to identify superdrop via integer
  //using IDType = EmptyID; // empty type of ID (for non-existent superdrop identity)

  size_t eps;    // multiplicity of superdroplet
  double radius; // radius of superdroplet
  double m_sol;  // mass of solute dissovled
  double coord3; // a 3rd spatial coordinate of superdroplet (z)
  double coord1; // a 1st spatial coordinate of superdroplet (x)
  double coord2; // a 2nd spatial coordinate of superdroplet (y)
  [[no_unique_address]] IDType id; // superdroplet (unique) identity

  Superdrop(const std::shared_ptr<const SoluteProperties> isolute,
            const size_t ieps, const double iradius,
            const double im_sol, const double icoord3,
            const double icoord1, const double icoord2,
            const IDType iid)
      : solute(isolute), eps(ieps), radius(iradius),
        m_sol(im_sol), coord3(icoord3), coord1(icoord1),
        coord2(icoord2), id(iid) {}

  inline double vol() const
  /* spherical volume of droplet calculated from its radius */
  {
    return 4.0 / 3.0 * M_PI * pow(radius, 3.0);
  }

  double mass() const;
  /* calculate total mass of droplet
  mass = (water + dry areosol)  */

  double superdroplet_wet_radius(const double s_ratio, const double temp) const;
  /* Performs Newton Raphson root finding algorithm using functions in 
  WetRadiusRootFinder struct to solve equation for equilibrium (wet) 
  radius of superdroplet at given relative humidity. Equilibrium radius 
  defined by radius when ODE from "An Introduction To Clouds...."
  (see note at top of file) eqn [7.28] = 0. */

  double akohler_factor(double temp) const;
  /* calculate value of a in raoult factor (exp^(a/r))
  to account for effect of dissolved solute
  on radial growth of droplet. Using equations from
  "An Introduction To Clouds...." (see note at top of file) */

  double bkohler_factor() const;
  /* calculate value of b in kelvin factor (1-b/r^3)
  to account for curvature on radial growth
  of droplet. Using equations from "An Introduction
  To Clouds...." (see note at top of file) */

  double change_radius(const double newradius);
  /* Update droplet radius to newradius or dry_radius() and
  return resultant change in radius. Prevents drops shrinking
  further once they are size of dry_radius(). */

  inline auto get_solute() const
  {
    return solute;
  }

  inline double get_dry_radius() const
  /* calculate radius as if dry droplet, ie.
  radius if drop was entirely made of solute */
  {
    return dry_radius();
  }
};

struct WetRadiusRootFinder
{
  private:
    double equilibrium_radius_polynomial(const double s_ratio, const double akoh,
                                        const double bkoh, const double ziter) const;
    /* returns value of (cubic) polynomial evaluted at ziter. Root of this
    polynomial is the value of the equilibrium (wet) radius of a superdroplet 
    at a given relative humidity (ie. s_ratio) derived from ODE in 
    "An Introduction To Clouds...." (see note at top of file) eqn [7.28] */                                    
    
    bool isnotconverged(const double new_ode, const double ode) const;
    /* boolean where True means criteria for ending newton raphson iterations
    has not yet been met. Criteria is standard local error test:
    |iteration - previous iteration| < RTOL * |iteration| + ATOL */
 
  public: 
    struct IterationReturn
    /* struct used for returning boolean and double
    from iterate_rootfinding_algorithm function */
    {
      bool do_iter;
      double ziter;
    };

    IterationReturn iterate_wetradius_rootfinding_algorithm(double ziter, const double s_ratio,
                                                            const double akoh, const double bkoh) const;
    /* performs 1 iteration of newton raphson root finding algorithm for 
    obtaining the equilibrium wet radius of the condensation ODE at a given 
    relative humidity (s_ratio). ODE from "An Introduction To Clouds...." 
    (see note at top of file) eqn [7.28] */    

};

struct SuperdropWithGbxindex
/* Structure containing a superdroplet (SD) and the index/unique
identifier of the gridbox (GBx) it occupies ie. the identity of
the GBx the SD is associated with */
{
  unsigned int sd_gbxindex; // index/unique identifier of gridbox the superdrop occupies
  Superdrop superdrop;

  SuperdropWithGbxindex(const unsigned int isd_gbxindex, Superdrop isuperdrop) 
      : sd_gbxindex(isd_gbxindex), superdrop(isuperdrop) {}
};

inline void sort_superdrops_via_gridboxindex(std::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* uses the value of sd_gbxindex within each SuperdropWithGbxindex
struct to sort the vector from lowest sd_gbxindex to highest. Sorting
of objects with same value of sd_gbxindex can take any order */
{
  auto compare = [](SuperdropWithGbxindex &a, SuperdropWithGbxindex &b)
  {
    return (a.sd_gbxindex) < (b.sd_gbxindex);
  };

  std::sort(SDsInGBxs.begin(), SDsInGBxs.end(), compare);
}

#endif // SUPERDROP_HPP