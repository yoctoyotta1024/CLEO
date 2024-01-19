/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: nullsuperdrops.hpp
 * Project: superdrops
 * Created Date: Friday 17th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 17th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions for handling null / empty
 * superdrops ie. with multiplicity, xi, = 0
 */

#ifndef NULLSUPERDROPS_HPP
#define NULLSUPERDROPS_HPP

#include <cassert>

#include <Kokkos_Core.hpp>
#include <Kokkos_NestedSort.hpp>

#include "../cleoconstants.hpp"
#include "./kokkosaliases_sd.hpp"
#include "./superdrop.hpp"

KOKKOS_INLINE_FUNCTION bool
is_null_superdrop(const Superdrop &drop)
/* raise error if multiplicity of
drop = 0, ie. superdrop is null */
{
  assert((drop.get_xi() > 0) && "superdrop xi < 1, null drop in coalescence");
  return 0;
}

KOKKOS_INLINE_FUNCTION subviewd_supers
is_null_supers(const subviewd_supers supers,
               const size_t nnull)
/* assert no null superdrops and return unchanged supers */
{
  assert((nnull == 0) && "no null superdrops should exist");
  return supers;
}

#endif // NULLSUPERDROPS_HPP
