/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: sortsupers.hpp
 * Project: gridboxes
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
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

namespace KCS = KokkosCleoSettings;

/* returns true if superdrops in "supers" view are already sorted by the Comparator */
template <class Comparator>
inline bool is_sorted_supers(const viewd_constsupers supers, const Comparator &comp) {
  return Kokkos::Experimental::is_sorted("IsSupersSorted", ExecSpace(), supers, comp);
}

/* Returns position of superdroplet count in counts/cumlcounts views given its
sdgbxindex. For in-domain superdroplets (0 <=sdgbxindex <= gbxindex_max),
position in counts/cumlcounts views is at the value of sdgbxindex, e.g.
if sdgbxindex=4, then position=4. If superdroplet has sdgbxindex > gbxindex_max
its position is last position of counts/cumlcounts array, i.e. all superdroplets with
sdgbxindex > gbxindex_max are found at last=(counts.extent(0) - 1) position of
counts/cumlcounts array. */
KOKKOS_INLINE_FUNCTION
size_t _get_count_position(const size_t sdgbxindex, const size_t gbxindex_max,
                           const viewd_counts counts) {
  return !(sdgbxindex <= gbxindex_max) ? (counts.extent(0) - 1) : sdgbxindex;
}

/* Functor used in parallel region of SortSupersBySdgbxindex "create_cumlcounts" (see below) */
struct CreateCumlcountsFunctor {
  size_t gbxindex_max;
  viewd_constsupers totsupers;
  viewd_counts counts;
  scatterviewd_counts s_counts;

  /* loop over superdroplets with this functor counts how many superdroplets have each gbxindex
  and how many have sdgbxindex > gbxindex_max. */
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk) const {
    auto counts_ = s_counts.access();
    const auto pos = _get_count_position(totsupers(kk).get_sdgbxindex(), gbxindex_max, counts);
    ++counts_(pos);  // atomic/duplicate "add" at ++counts(sdgbxindex) or ++counts([last]);
  };
};

/* Functor used in parallel regions of SortSupersBySdgbxindex counting_sort functions (see below) */
struct CountingSortFunctor {
  size_t gbxindex_max;
  subviewd_supers supers;
  viewd_supers totsupers_tmp;
  viewd_counts cumlcounts;

  /* loop over superdroplets with this functor copies the superdroplets from supers to totsupers_tmp
  in new order based on their sdgbxindexes such that they are sorted from lowest to highest
  gbxindex (for superdroplets with sdgbxindex <= gbxindex_max). In the sorted view, superdroplets
  with sdgbxindex > gbxindex_max occur after those with sdgbxindex <= gbxindex_max but may not be
  sorted amongst themselves. sdgbxindex is set to LIMITVALUES::oob_gbxindex for all superdroplets in
  "supers" that are copied/moved to totsupers_tmp. */
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk) const {
    const auto pos = _get_count_position(supers(kk).get_sdgbxindex(), gbxindex_max, cumlcounts);
    auto new_kk = Kokkos::atomic_fetch_add(&cumlcounts(pos), 1);
    totsupers_tmp(new_kk) = supers(kk);
    supers(kk).set_sdgbxindex(LIMITVALUES::oob_gbxindex);  // fail-safe reset totsupers
  }
};

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

  /* counts number of superdroplets in each gbx with sdgbxindex <= gbxindex_max and all
  superdroplets with sdgbxindex > gbxindex_max. E.g. if totsupers contains 5 supersdroplets
  with sdgbxindex=0, then counts array at position 0 will be 5, meanwhile all counts of
  superdroplets with sdgbxindex > gbxindex_max go into last position of counts array.
  Returns cumulative sum from position 0 to last position.
  TODO(CB): (optional) improvement when scatterview uses atomics (e.g. on GPUs) would be to loop
  over domainsupers and oobsupers like couting_sort(...) can. */
  viewd_counts create_cumlcounts(const viewd_constsupers totsupers) {
    const auto ntotsupers = size_t{totsupers.extent(0)};
    s_counts.reset();
    const auto functor = CreateCumlcountsFunctor{gbxindex_max, totsupers, counts, s_counts};
    Kokkos::parallel_for("increment_counts", Kokkos::RangePolicy<ExecSpace>(0, ntotsupers),
                         functor);
    Kokkos::Experimental::contribute(counts, s_counts);

    Kokkos::Experimental::exclusive_scan("cumulative_sum", ExecSpace(), counts, cumlcounts, 0);
    Kokkos::Experimental::fill(ExecSpace(), counts, 0);  // reset in preparation for next call

    return cumlcounts;
  }

  /* Part of counting sort algorithm involving copying/movement of superdroplets into new sorted
  array. Takes superdroplets from start to end of original totsupers view and copies them to
  their positions in a new view according to the CountingSortFunctor. May also modify superdroplets
  in totsupers, e.g. setting their sdgbxindex to LIMITVALUES::oob_gbxindex. */
  viewd_supers counting_sort(const viewd_supers totsupers) {
    const auto ntotsupers = size_t{totsupers.extent(0)};
    const auto functor = CountingSortFunctor{gbxindex_max, totsupers, totsupers_tmp, cumlcounts};
    Kokkos::parallel_for("counting_sort", Kokkos::RangePolicy<ExecSpace>(0, ntotsupers), functor);

    /* assertion for debugging only works for hostspace cumlcounts */
    // assert((cumlcounts(counts.extent(0) - 1) == totsupers_tmp.extent(0)) &&
    //        "last cumulative sum of totsupers count should equal expected number of totsupers");

    return totsupers_tmp;
  }

  /* same output as counting_sort(totsupers) but uses heirarchal loop over d_gbxs to reduce chance
  of in atomic operation. In doing so, it implicitly totsupers=domainsupers+oob_supers,
  i.e. that domainsupers is a subview of totsupers which starts and the same address as
  totsupers and that oob_supers is a subview which starts at the end of domainsupers and ends
  at the end of totsupers. */
  viewd_supers counting_sort(const viewd_constgbx d_gbxs, const subviewd_supers domainsupers,
                             const subviewd_supers oob_supers) {
    const auto ngbxs = size_t{d_gbxs.extent(0)};
    Kokkos::parallel_for(
        "counting_sort_gbxs", TeamPolicy(ngbxs, KCS::team_size),
        KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
          const auto ii = team_member.league_rank();
          auto supers = d_gbxs(ii).supersingbx(domainsupers);
          const auto nsupers = size_t{supers.extent(0)};
          const auto functor = CountingSortFunctor{gbxindex_max, supers, totsupers_tmp, cumlcounts};
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, nsupers), functor);
        });

    const auto noobs = size_t{oob_supers.extent(0)};
    const auto functor = CountingSortFunctor{gbxindex_max, oob_supers, totsupers_tmp, cumlcounts};
    Kokkos::parallel_for("counting_sort_oob", Kokkos::RangePolicy<ExecSpace>(0, noobs), functor);

    /* assertion for debugging only works for hostspace cumlcounts */
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

  /* Counting sort algorithm to (stable) sort superdroplets inside the domain by sdgbxindex.
  Superdrops in totsupers may change (e.g. sdgbxindex may be set to LIMITVALUES::oob_gbxindex) and
  returned view may not be the same totsupers view given as argument. Superdroplets outside of the
  domain (i.e. sdgbxindex > gbxindex_max) are not guarenteed to be sorted. This operator calls
  sorting functions with loop over gridboxes in order to reduce the chance of conflicts in atomic
  operations. In doing so it assumes totsupers=domainsupers+oob_supers, i.e. that domainsupers is
  a subview of totsupers which starts and the same address as totsupers and that oob_supers is a
  subview which starts at the end of domainsupers and ends at the end of totsupers.
  */
  viewd_supers operator()(const viewd_supers totsupers, const viewd_constgbx d_gbxs,
                          kkpair_size_t domainrefs) {
    const auto domainsupers = Kokkos::subview(totsupers, domainrefs);
    const kkpair_size_t oobrefs = Kokkos::make_pair(domainrefs.second, totsupers.extent(0));
    const auto oob_supers = Kokkos::subview(totsupers, oobrefs);

    cumlcounts = create_cumlcounts(totsupers);
    const auto sorted_supers = counting_sort(d_gbxs, domainsupers, oob_supers);
    totsupers_tmp = totsupers;  // fail-safe reset totsupers_tmp
    return sorted_supers;
  }
};

#endif  // LIBS_GRIDBOXES_SORTSUPERS_HPP_
