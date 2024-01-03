/*
 * ----- CLEO -----
 * File: coalbure_flag.hpp
 * Project: superdrops
 * Created Date: Friday 29th December 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 3rd January 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * concept and structs that return a flag used
 * in DoCoalBuRe to decide which out of breakup,
 * coalescence or rebound should occur.
 */

#ifndef COALBURE_FLAG_HPP
#define COALBURE_FLAG_HPP

#include <functional>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "./collisionkinetics.hpp"
#include "./superdrop.hpp"
#include "./terminalvelocity.hpp"

template <typename F>
concept CoalBuReFlag = requires(F f,
                                const Superdrop &d1,
                                const Superdrop &d2)
/* operator returns flag indicating rebound or
coalescence or breakup. If flag = 1 -> coalescence.
If flag = 2 -> breakup. Otherwise -> rebound. */
{
  {
    f(d1, d2)
  } -> std::convertible_to<unsigned int>;
};

struct SUCoalBuReFlag
{
  KOKKOS_FUNCTION
  unsigned int operator()(const Superdrop &drop1,
                          const Superdrop &drop2) const;
  /* function returns flag indicating rebound or
  coalescence or breakup. If flag = 1 -> coalescence.
  If flag = 2 -> breakup. Otherwise -> rebound.
  Flag decided based on the kinetic arguments in
  section 2.2 of Szakáll and Urbich 2018
  (neglecting grazing angle considerations) */
};

struct TSCoalBuReFlag
{
  KOKKOS_FUNCTION
  unsigned int operator()(const Superdrop &drop1,
                          const Superdrop &drop2) const;
  /* function returns flag indicating rebound or
  coalescence or breakup. If flag = 1 -> coalescence.
  If flag = 2 -> breakup. Otherwise -> rebound.
  Flag decided based on the kinetic arguments from
  section 4 of Testik et al. 2011 (figure 12) as well 
  as coalescence efficiency from Straub et al. 2010 */
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION unsigned int
SUCoalBuReFlag::operator()(const Superdrop &drop1,
                              const Superdrop &drop2) const
/*  function returns flag indicating rebound or
coalescence or breakup. If flag = 1 -> coalescence.
If flag = 2 -> breakup. Otherwise -> rebound.
Flag decided based on the kinetic arguments in
section 2.2 of Szakáll and Urbich 2018
(neglecting grazing angle considerations) */
{
  const auto r1 = drop1.get_radius();
  const auto r2 = drop2.get_radius();
  const auto terminalv = SimmelTerminalVelocity{};

  const auto cke = collision_kinetic_energy(r1, r2,
                                            terminalv(drop1),
                                            terminalv(drop2));

  if (cke < surfenergy(Kokkos::fmin(r1, r2))) // cke < surface energy of small drop
  {
    return 0; // rebound
  }
  else if (cke < coal_surfenergy(r1, r2)) // Weber number < 1
  {
    return 1; // coalescence
  }
  else // Weber number > 1
  {
    return 2; // breakup
  }
}

KOKKOS_FUNCTION unsigned int
TSCoalBuReFlag::operator()(const Superdrop &drop1,
                              const Superdrop &drop2) const
/* function returns flag indicating rebound or
coalescence or breakup. If flag = 1 -> coalescence.
If flag = 2 -> breakup. Otherwise -> rebound.
Flag decided based on the kinetic arguments from
section 4 of Testik et al. 2011 (figure 12) as well
as coalescence efficiency from Straub et al. 2010 */
{
  const auto r1 = drop1.get_radius();
  const auto r2 = drop2.get_radius();
  const auto terminalv = SimmelTerminalVelocity{};

  const auto cke = collision_kinetic_energy(r1, r2,
                                            terminalv(drop1),
                                            terminalv(drop2));

  if (cke < surfenergy(Kokkos::fmin(r1, r2))) // cke < surface energy of small drop
  {
    return coalescence_or_rebound(); // below DE2 boundary
  }
  else if (cke < surfenergy(Kokkos::fmax(r1, r2))) // cke < surface energy of large drop
  {
    return coalescence_or_breakup(); // below DE1 boundary
  }
  else // above DE1 boundary
  {
    return 2; // breakup
  }
}

#endif // COALBURE_FLAG_HPP
