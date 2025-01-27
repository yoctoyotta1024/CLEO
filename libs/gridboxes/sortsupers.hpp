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
// #include <cassert>

#include "../kokkosaliases.hpp"
#include "superdrops/superdrop.hpp"

/* returns true if superdrops in "supers" view are already sorted by the Comparator */
template <class Comparator>
inline bool is_sorted_supers(const viewd_constsupers supers, const Comparator &comp) {
  return Kokkos::Experimental::is_sorted("IsSupersSorted", ExecSpace(), supers, comp);
}

/* Counting sort algorithm to (stable) sort superdroplets inside the domain by sdgbxindex.
Gridbox indexes are assumed to run from 0 to gbxindex_max so that Superdroplets inside the
domain have 0 <= sdgbxindex <= gbxindex_max. Superdroplets outside of the domain
(i.e. sdgbxindex > gbxindex_max) are not guarenteed to be sorted. */
struct SortSupersBySdgbxindex {
 public:
  size_t gbxindex_max;          /**< maximum gbxindex of in-domain superdroplets */
  viewd_counts counts;          /**< number of superdroplets in each gridbox + outside of domain */
  scatterviewd_counts s_counts; /**< scatter view for atomics/duplicate ops with counts */
  viewd_counts cumlcounts;      /**< cumulative version of counts */
  viewd_supers totsupers_tmp;   /**< temporary view of superdroplets used by sorting algorithm */

  SortSupersBySdgbxindex(const size_t gbxindex_max_, const size_t ntotsupers)
      : gbxindex_max(gbxindex_max_),
        counts("counts", gbxindex_max + 1),
        s_counts(counts),
        cumlcounts("cumlcounts", gbxindex_max + 1),
        totsupers_tmp("totsupers_tmp", ntotsupers) {}

  /* a precedes b if its sdgbxindex is smaller */
  struct SortComparator {
    KOKKOS_INLINE_FUNCTION
    bool operator()(const Superdrop &a, const Superdrop &b) const {
      return (a.get_sdgbxindex()) < (b.get_sdgbxindex());
    }
  };

  /* returns true if superdrops in view are sorted by their sdgbxindexes in ascending order */
  bool is_sorted(const viewd_constsupers supers) const {
    return is_sorted_supers(supers, SortComparator{});
  }

  /* Returns position of superdroplet count in counts/cumlcounts views given its
  sdgbxindex. For in-domain superdroplets (0 <=sdgbxindex <= gbxindex_max),
  position in counts/cumlcounts views is at the value of sdgbxindex, e.g.
  if sdgbxindex=4, then position=4. If superdroplet has sdgbxindex > gbxindex_max
  its position is last position of counts/cumlcounts array, i.e. all superdroplets with
  sdgbxindex > gbxindex_max are found at last=(counts.extent(0) - 1) position of
  counts/cumlcounts array. */
  KOKKOS_INLINE_FUNCTION
  size_t get_count_position(const size_t sdgbxindex) const {
    return !(sdgbxindex <= gbxindex_max) ? (counts.extent(0) - 1) : sdgbxindex;
  }

  /* counts number of superdroplets in each gbx with sdgbxindex <= gbxindex_max and all
  superdroplets with sdgbxindex > gbxindex_max. E.g. if totsupers contains 5 supersdroplets
  with sdgbxindex=0, then counts array at position 0 will be 5, meanwhile all counts of
  superdroplets with sdgbxindex > gbxindex_max go into last position of counts array.
  Returns cumulative sum from position 0 to last position.  */
  viewd_counts create_cumlcounts(const viewd_constsupers totsupers) {
    const auto ntotsupers = size_t{totsupers.extent(0)};
    s_counts.reset();
    Kokkos::parallel_for(
        "increment_counts", Kokkos::RangePolicy<ExecSpace>(0, ntotsupers),
        KOKKOS_CLASS_LAMBDA(const size_t kk) {
          auto counts_ = s_counts.access();
          auto pos = get_count_position(totsupers(kk).get_sdgbxindex());
          ++counts_(pos);  // atomic/duplicate "add" at ++counts(sdgbxindex) or ++counts([last]);
        });

    Kokkos::Experimental::exclusive_scan("cumulative_sum", ExecSpace(), counts, cumlcounts, 0);
    Kokkos::Experimental::fill(ExecSpace(), counts, 0);  // reset in preparation for next call

    return cumlcounts;
  }

  viewd_supers counting_sort(const viewd_supers totsupers) {
    const auto ntotsupers = size_t{totsupers.extent(0)};
    Kokkos::parallel_for(
        "counting_sort", Kokkos::RangePolicy<ExecSpace>(0, ntotsupers),
        KOKKOS_CLASS_LAMBDA(const size_t kk) {
          const auto pos = get_count_position(totsupers(kk).get_sdgbxindex());
          auto new_kk = Kokkos::atomic_fetch_add(&cumlcounts(pos), 1);
          totsupers_tmp(new_kk) = totsupers(kk);
          totsupers(kk).set_sdgbxindex(LIMITVALUES::oob_gbxindex);  // fail-safe reset totsupers
        });

    // /* assertion for debugging only works for hostspace cumlcounts */
    // assert((cumlcounts(counts.extent(0) - 1) == totsupers_tmp.extent(0)) &&
    //        "last cumulative sum of totsupers count should equal expected number of totsupers");

    return totsupers_tmp;
  }

  /* Counting sort algorithm to (stable) sort superdroplets inside the domain by sdgbxindex.
  Superdrops in totsupers may change (e.g. sdgbxindex may be set to LIMITVALUES::oob_gbxindex) and
  returned view may not be the same totsupers view given as argument. Superdroplets outside of the
  domain (i.e. sdgbxindex > gbxindex_max) are not guarenteed to be sorted. */
  viewd_supers operator()(const viewd_supers totsupers) {
    cumlcounts = create_cumlcounts(totsupers);
    const auto sorted_supers = counting_sort(totsupers);
    totsupers_tmp = totsupers;  // fail-safe reset totsupers_tmp
    return sorted_supers;
  }
};

#endif  // LIBS_GRIDBOXES_SORTSUPERS_HPP_
