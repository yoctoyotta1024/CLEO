/*
 * ----- CLEO -----
 * File: sortsupers.hpp
 * Project: gridboxes
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 19th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions used when sorting /shuffling superdrops
 * e.g. based on their gridbox indexes
 */

#ifndef SORTSUPERS_HPP
#define SORTSUPERS_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_NestedSort.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "../kokkosaliases.hpp"
#include "superdrops/superdrop.hpp"

struct SortComparator
/* a precedes b if its sdgbxindex is smaller */
{
  KOKKOS_INLINE_FUNCTION
  bool operator()(const Superdrop &a, const Superdrop &b) const
  {
    return (a.get_sdgbxindex()) < (b.get_sdgbxindex());
  }
};

inline viewd_supers sort_supers(const viewd_supers supers)
/* sort a view of superdroplets by their sdgbxindexes
so that superdrops in the view are ordered from
lowest to highest sdgbxindex. Note that sorting of
superdrops with matching sdgbxindex can take any order */
{
  using TeamPol = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;

  Kokkos::parallel_for(
      "sortingsupers", TeamPol(1, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const TeamPol::member_type &t) {
        Kokkos::Experimental::sort_team(t, supers, SortComparator{});
      });

  return supers;
}

inline bool is_sorted(const viewd_constsupers supers)
/* returns true if superdrops in supers view are
sorted by their sdgbxindexes in ascending order */
{
  return Kokkos::Experimental::
      is_sorted("IsSupersSorted",
                Kokkos::DefaultExecutionSpace(),
                supers,
                SortComparator{});
}

#endif // SORTSUPERS_HPP
