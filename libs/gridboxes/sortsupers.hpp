/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: sortsupers.hpp
 * Project: gridboxes
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions used when sorting /shuffling superdrops
 * e.g. based on their gridbox indexes
 */

#ifndef LIBS_GRIDBOXES_SORTSUPERS_HPP_
#define LIBS_GRIDBOXES_SORTSUPERS_HPP_

#include <algorithm>

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "../kokkosaliases.hpp"
#include "superdrops/superdrop.hpp"

/* a precedes b if its sdgbxindex is smaller */
struct SortComparator {
  KOKKOS_INLINE_FUNCTION
  bool operator()(const Superdrop &a, const Superdrop &b) const {
    return (a.get_sdgbxindex()) < (b.get_sdgbxindex());
  }
};

/* sort a view of superdroplets by their sdgbxindexes
so that superdrops in the view are ordered from
lowest to highest sdgbxindex. Note that sorting of
superdrops with matching sdgbxindex can take any order */
inline viewd_supers sort_supers(const viewd_supers supers) {
  Kokkos::sort(ExecSpace(), supers, SortComparator{});

  return supers;
}

/* returns true if superdrops in supers view are
sorted by their sdgbxindexes in ascending order */
inline bool is_sorted(const viewd_constsupers supers) {
  return Kokkos::Experimental::is_sorted("IsSupersSorted", ExecSpace(), supers, SortComparator{});
}

#endif  // LIBS_GRIDBOXES_SORTSUPERS_HPP_
