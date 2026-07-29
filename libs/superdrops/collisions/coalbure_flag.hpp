/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: coalbure_flag.hpp
 * Project: collisions
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept and structs that return a flag used
 * in DoCoalBuRe to decide whether breakup,
 * coalescence or rebound should occur.
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_COALBURE_FLAG_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_COALBURE_FLAG_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <functional>
#include <random>

#include "../superdrop.hpp"
#include "../terminalvelocity.hpp"
#include "./coalescence_efficiency.hpp"
#include "./collisionkinetics.hpp"

/* operator returns flag indicating rebound or
coalescence or breakup. If flag = 1 -> coalescence.
If flag = 2 -> breakup. Otherwise -> rebound. */
template <typename F>
concept CoalBuReFlag = requires(F f, const double phi, const Superdrop& d1, const Superdrop& d2) {
  { f(phi, d1, d2) } -> std::convertible_to<unsigned int>;
};

template <VelocityFormula TerminalVelocity>
struct SUCoalBuReFlag {
 private:
  TerminalVelocity terminalv;

  /* function returns flag indicating rebound or coalescence or breakup.
  If flag = 1 -> coalescence.
  If flag = 2 -> breakup.
  Otherwise -> rebound.
  Flag decided based on the kinetic arguments in section 2.2 of
  Szakáll and Urbich 2018 (neglecting grazing angle considerations) */
  KOKKOS_FUNCTION
  unsigned int operator()(const Superdrop& drop1, const Superdrop& drop2) const {
    const auto r1 = drop1.get_radius();
    const auto r2 = drop2.get_radius();

    const auto cke = collision_kinetic_energy(r1, r2, terminalv(drop1),
                                              terminalv(drop2));  // [J]

    if (cke < surfenergy(Kokkos::fmin(r1, r2))) {  // cke < surface energy of small drop
      return 0;                                    // rebound
    } else if (cke < coal_surfenergy(r1, r2)) {    // Weber number < 1
      return 1;                                    // coalescence
    } else {                                       // Weber number > 1
      return 2;                                    // breakup
    }
  }

 public:
  explicit SUCoalBuReFlag(TerminalVelocity tv) : terminalv(tv) {}

  /* adaptor of operator to satisfy CoalBuReFlag concept */
  KOKKOS_FUNCTION
  unsigned int operator()(const double phi, const Superdrop& drop1, const Superdrop& drop2) const {
    return operator()(drop1, drop2);
  }
};

template <VelocityFormula TerminalVelocity>
struct TSCoalBuReFlag {
 private:
  TerminalVelocity terminalv;

  /* returns flag that indicates coalescence (flag=1) or rebound (flag=0)
  based on coalescence efficiency from Straub et al. 2010 */
  KOKKOS_FUNCTION
  unsigned int rebound_or_coalescence(const Superdrop& drop1, const Superdrop& drop2,
                                      const double phi, const double cke) const {
    if (is_coalescence(drop1, drop2, phi, cke)) {
      return 1;  // coalescence
    } else {
      return 0;  // rebound
    }
  }

  /* returns flag that indicates coalescence (flag=1) or breakup (flag=2)
  based on coalescence efficiency from Straub et al. 2010 */
  KOKKOS_FUNCTION
  unsigned int coalescence_or_breakup(const Superdrop& drop1, const Superdrop& drop2,
                                      const double phi, const double cke) const {
    if (is_coalescence(drop1, drop2, phi, cke)) {
      return 1;  // coalescence
    } else {
      return 2;  // breakup
    }
  }

  /* returns true if comparison of random number with coalescence efficiency from
  Straub et al. 2010 indicates coalescence should occur */
  KOKKOS_FUNCTION bool is_coalescence(const Superdrop& drop1, const Superdrop& drop2,
                                      const double phi, const double cke) const {
    const auto ecoal = coalescence_efficiency_straub2010(drop1, drop2, cke);

    if (phi < ecoal) {
      return true;
    } else {
      return false;
    }
  }

 public:
  explicit TSCoalBuReFlag(TerminalVelocity tv) : terminalv(tv) {}

  /* function returns flag indicating rebound or coalescence or breakup.
  If flag = 1 -> coalescence.
  If flag = 2 -> breakup.
  Otherwise -> rebound.
  Flag decided based on the kinetic arguments from section 4 of Testik et al. 2011 (figure 12;
  first proposed in Testik 2009) as well as the coalescence efficiency from Straub et al. 2010 */
  KOKKOS_FUNCTION
  unsigned int operator()(const double phi, const Superdrop& drop1, const Superdrop& drop2) const {
    const auto r1 = drop1.get_radius();
    const auto r2 = drop2.get_radius();

    const auto cke = collision_kinetic_energy(r1, r2, terminalv(drop1), terminalv(drop2));

    if (cke < surfenergy(Kokkos::fmin(r1, r2))) {             // cke < surface energy of small drop
      return rebound_or_coalescence(drop1, drop2, phi, cke);  // below DE2 boundary
    } else if (cke < surfenergy(Kokkos::fmax(r1, r2))) {      // cke < surface energy of large drop
      return coalescence_or_breakup(drop1, drop2, phi, cke);  // below DE1 boundary
    } else {                                                  // above DE1 boundary
      return 2;                                               // breakup
    }
  }
};

template <VelocityFormula TerminalVelocity>
struct StraubCoalBuReFlag {
 private:
  TerminalVelocity terminalv;

 public:
  explicit StraubCoalBuReFlag(TerminalVelocity tv) : terminalv(tv) {}

  /* function returns flag indicating coalescence or breakup (never rebound).
   * If flag = 1 -> coalescence. If flag = 2 -> breakup.
   * Flag decided based on Straub et al. 2010 coalescence efficiency equation 3
   * (with definition of terms from Schlottle et al. 2010) compared to random number "phi"
   * function signature matches conditions to satisfy CoalBuReFlag concept
   * */
  KOKKOS_INLINE_FUNCTION
  unsigned int operator()(const double phi, const Superdrop& drop1, const Superdrop& drop2) const {
    const auto r1 = drop1.get_radius();
    const auto r2 = drop2.get_radius();
    const auto cke = collision_kinetic_energy(r1, r2, terminalv(drop1),
                                              terminalv(drop2));  // [J]

    const auto ecoal = coalescence_efficiency_straub2010(drop1, drop2, cke);

    if (phi < ecoal) {
      return 1;  // coalescence
    } else {
      return 2;  // breakup
    }
  }
};

struct ConstCoalBuReFlag {
 private:
  double coaleff;  // flag indicating whether coalescence, breakup or rebound should occur
 public:
  explicit ConstCoalBuReFlag(const double coaleff) : coaleff(coaleff) {}

  /* function returns flag indicating coalescence or breakup (never rebound).
   * If flag = 1 -> coalescence. If flag = 2 -> breakup.
   * Flag decided based on constant coalescence efficiency compared to random number "phi"
   * function signature matches conditions to satisfy CoalBuReFlag concept
   * */
  KOKKOS_INLINE_FUNCTION
  unsigned int operator()(const double phi, const Superdrop& drop1, const Superdrop& drop2) const {
    if (phi < coaleff) {
      return 1;  // coalescence
    } else {
      return 2;  // breakup
    }
  }
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_COALBURE_FLAG_HPP_
