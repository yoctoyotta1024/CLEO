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
#include "./collisionkinetics.hpp"

/* operator returns flag indicating rebound or
coalescence or breakup. If flag = 1 -> coalescence.
If flag = 2 -> breakup. Otherwise -> rebound. */
template <typename F>
concept CoalBuReFlag = requires(F f, const double phi, const Superdrop &d1, const Superdrop &d2) {
  { f(phi, d1, d2) } -> std::convertible_to<unsigned int>;
};

struct SUCoalBuReFlag {
 private:
  /* function returns flag indicating rebound or
  coalescence or breakup. If flag = 1 -> coalescence.
  If flag = 2 -> breakup. Otherwise -> rebound.
  Flag decided based on the kinetic arguments in
  section 2.2 of Szakáll and Urbich 2018
  (neglecting grazing angle considerations) */
  KOKKOS_FUNCTION
  unsigned int operator()(const Superdrop &drop1, const Superdrop &drop2) const;

 public:
  /* adaptor of operator to satisfy CoalBuReFlag concept */
  KOKKOS_FUNCTION
  unsigned int operator()(const double phi, const Superdrop &drop1, const Superdrop &drop2) const {
    return operator()(drop1, drop2);
  }
};

struct TSCoalBuReFlag {
 private:
  /* returns flag that indicates coalescence (flag=1)
  or rebound (flag=0) based on coalescence efficiency
  from Straub et al. 2010 */
  KOKKOS_FUNCTION
  unsigned int rebound_or_coalescence(const Superdrop &drop1, const Superdrop &drop2,
                                      const double phi, const double cke) const;

  /* returns flag that indicates coalescence (flag=1)
  or breakup (flag=2) based on coalescence efficiency
  from Straub et al. 2010 */
  KOKKOS_FUNCTION
  unsigned int coalescence_or_breakup(const Superdrop &drop1, const Superdrop &drop2,
                                      const double phi, const double cke) const;

  /* returns truw if comparison of random numnber
  with coalescence efficiency from Straub et al. 2010
  indicates coalescence should occur */
  KOKKOS_FUNCTION bool is_coalescence(const Superdrop &drop1, const Superdrop &drop2,
                                      const double phi, const double cke) const;

  /* coalescence efficency given a collision occurs
  according to parameterisation from Straub et al. 2010
  section 3, equation 5 and Schlottke et al. 2010
  section 4a equation 11 */
  KOKKOS_FUNCTION
  double coalescence_efficiency(const Superdrop &drop1, const Superdrop &drop2,
                                const double cke) const;

 public:
  TSCoalBuReFlag() {}

  /* function returns flag indicating rebound or
  coalescence or breakup. If flag = 1 -> coalescence.
  If flag = 2 -> breakup. Otherwise -> rebound.
  Flag decided based on the kinetic arguments from
  section 4 of Testik et al. 2011 (figure 12) as well
  as coalescence efficiency from Straub et al. 2010 */
  KOKKOS_FUNCTION
  unsigned int operator()(const double phi, const Superdrop &drop1, const Superdrop &drop2) const;
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_COALBURE_FLAG_HPP_
