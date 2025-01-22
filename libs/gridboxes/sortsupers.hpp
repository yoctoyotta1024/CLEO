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

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "../kokkosaliases.hpp"
#include "superdrops/superdrop.hpp"

/* sort the "supers" view of superdroplets by the Comparator. Note that sorting of superdrops
is not guarenteed to be stable (not guarenteed to maintain previous order of equal values). */
template <class Comparator>
inline viewd_supers sort_supers(const viewd_supers supers, const Comparator &comp) {
  Kokkos::sort(ExecSpace(), supers, comp);
  return supers;
}

/* returns true if superdrops in "supers" view are already sorted by the Comparator */
template <class Comparator>
inline bool is_sorted_supers(const viewd_constsupers supers, const Comparator &comp) {
  return Kokkos::Experimental::is_sorted("IsSupersSorted", ExecSpace(), supers, comp);
}

struct SortSupersBySdgbxindex {
 private:
  /* a precedes b if its sdgbxindex is smaller */
  struct SortComparator {
    KOKKOS_INLINE_FUNCTION
    bool operator()(const Superdrop &a, const Superdrop &b) const {
      return (a.get_sdgbxindex()) < (b.get_sdgbxindex());
    }
  };

 public:
  SortSupersBySdgbxindex() {}

  /* sort the "supers" view of superdroplets by their sdgbxindexes so that superdrops in the
  view are ordered from lowest to highest sdgbxindex. Note that sorting of superdrops with matching
  sdgbxindex is not guarenteed to maintain their previous order (unstable sorting) */
  viewd_supers operator()(const viewd_supers supers) {
    return sort_supers(supers, SortComparator{});
  }

  /* returns true if superdrops in view are sorted by their sdgbxindexes in ascending order */
  bool is_sorted(const viewd_constsupers supers) const {
    return is_sorted_supers(supers, SortComparator{});
  }
};

#endif  // LIBS_GRIDBOXES_SORTSUPERS_HPP_
