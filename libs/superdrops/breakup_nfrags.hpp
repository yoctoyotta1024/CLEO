/*
 * ----- CLEO -----
 * File: breakup_nfrags.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 23rd November 2023 
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * concept and structures for calculating
 * the number of fragments produce.
 * Used e.g. by the DoBreakup struct.
*/

#ifndef BREAKUP_NFRAGS_HPP
#define BREAKUP_NFRAGS_HPP

#include <concepts>

#include <Kokkos_Core.hpp>

template <typename F>
concept NFragments = requires(F f,
                              const Superdrop &d1,
                              const Superdrop &d2)
/* Objects that are of type 'NFragments'
take a pair of superdroplets and returns
something convertible to a double (such as
the number of fragments from a breakup event) */
{
  {
    f(d1, d2)
  } -> std::convertible_to<double>;
};

struct ConstNFrags
/* operator always returns constant 
number of fragments 'nfrags'. Struct
obeys NFragments concept  */
{
private:
  const double nfrags;
public:
  ConstNFrags(const double nfrags) : nfrags(nfrags) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &d1,
                    const Superdrop &d2) const
  /* always returns constant number of fragments 'nfrags' */
  {
    return nfrags;
  }
};

#endif // BREAKUP_NFRAGS_HPP