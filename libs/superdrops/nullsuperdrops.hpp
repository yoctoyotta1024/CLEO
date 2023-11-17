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

// KOKKOS_INLINE_FUNCTION bool
// if_null_superdrop(Superdrop &drop)
// /* if multiplicity of drop = 0, ie. superdrop
// is null, so change it's sdgbxindex to be value that
// indicates superdrop is out of domain (ie. no
// longer exists) and return true */
// {
//   if (drop.get_xi() < 1) // ie. xi == 0
//   {
//     drop.set_sdgbxindex(LIMITVALUES::uintmax);
//     return 1;
//   }
//   else
//   {
//     return 0;
//   }
// }

// KOKKOS_INLINE_FUNCTION subviewd_supers
// remove_null_supers(subviewd_supers supers,
//                    const size_t nnull)
// /* sort view of superdroplets by their sdgbxindexes
// from lowest to highest sdgbxindex. Then set new
// subview excluding 'nnull' number of supers */
// {
//   using TeamPol = Kokkos::TeamPolicy<ExecSpace>;

//   struct SortComparator
//   /* a precedes b if its sdgbxindex is smaller */
//   {
//     KOKKOS_INLINE_FUNCTION
//     bool operator()(const Superdrop &a, const Superdrop &b) const
//     {
//       return (a.get_sdgbxindex()) < (b.get_sdgbxindex());
//     }
//   };

//   Kokkos::parallel_for(
//       "sortingsupers_thread", TeamPol(1, Kokkos::AUTO()),
//       KOKKOS_LAMBDA(const TeamPol::member_type &t) {
//         Kokkos::Experimental::sort_thread(t, supers, SortComparator{});
//       });

//   const size_t nsupers(supers.extent(0) - nnull); // exclude null supers from no. supers
//   const Kokkos::pair<size_t, size_t> new_refs({0, nsupers});

//   return Kokkos::subview(supers, new_refs);
// }

#endif // NULLSUPERDROPS_HPP 

