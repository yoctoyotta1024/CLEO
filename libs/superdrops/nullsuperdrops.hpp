/*
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
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions for handling null / empty 
 * superdrops ie. with multiplicity, xi, = 0
 */

#ifndef NULLSUPERDROPS_HPP 
#define NULLSUPERDROPS_HPP 

#include <cassert>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "./superdrop.hpp"

KOKKOS_INLINE_FUNCTION void
is_null_superdrop(const Superdrop &drop)
/* raise error if multiplicity of
drop = 0, ie. superdrop is null */
{
  assert((drop.get_xi() > 0) && "superdrop xi < 1, null drop in coalescence");
}

KOKKOS_INLINE_FUNCTION void
if_null_superdrop(Superdrop &drop)
/* if multiplicity of drop = 0, ie. superdrop
is null, raise error or (uncomment if desired)
change it's sdgbxindex to be value that indicates
superdrop is out of domain (ie. no longer exists) */
{
  if (drop.get_xi() < 1) // ie. xi == 0
  {
    drop.set_sdgbxindex(LIMITVALUES::uintmax);
  }
}

#endif // NULLSUPERDROPS_HPP 

