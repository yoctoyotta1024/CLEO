/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: supersindomain.hpp
 * Project: gridboxes
 * Created Date: Tuesday 21st January 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 21st January 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions and structures related to handling superdroplets inside CLEO's domain (on one node)
 */

#ifndef LIBS_GRIDBOXES_SUPERSINDOMAIN_HPP_
#define LIBS_GRIDBOXES_SUPERSINDOMAIN_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "gridboxes/findrefs.hpp"
#include "gridboxes/sortsupers.hpp"
#include "superdrops/kokkosaliases_sd.hpp"

/* References to identify the chunk of memory containing super-droplets occupying domain
(i.e. within any of the gridboxes on a single node), e.g. through std::span or Kokkos::subview) */
struct SupersInDomain {
 private:
  Kokkos::pair<unsigned int, unsigned int> gbxindex_range; /**< {min, max} gbxindex of domain */
  viewd_supers totsupers; /**< view of all superdrops (both in and out of bounds of domain) */
  kkpair_size_t
      domainrefs; /**< position in view of (first, last) superdrop that occupies gridbox */
  SortSupersBySdgbxindex sort_by_sdgbxindex; /**< method to sort view of superdrops by sdgbxindex */

  void set_totsupers_domainrefs(const viewd_supers totsupers_) {
    totsupers = totsupers_;
    domainrefs = find_domainrefs(totsupers, gbxindex_range);
  }

 public:
  SupersInDomain() = default;   // Kokkos requirement for a (dual)View
  ~SupersInDomain() = default;  // Kokkos requirement for a (dual)View

  explicit SupersInDomain(const viewd_supers totsupers_, const unsigned int gbxindex_max)
      : gbxindex_range({0, gbxindex_max}),
        totsupers(totsupers_),
        domainrefs({0, 0}),
        sort_by_sdgbxindex(SortSupersBySdgbxindex()) {
    totsupers = sort_by_sdgbxindex(totsupers_);
    set_totsupers_domainrefs(totsupers_);
  }

  viewd_supers get_totsupers() const {
    return totsupers;
  }  // TODO(CB): replace with appending SDs func which could resize(?)

  viewd_constsupers get_totsupers_readonly() const { return totsupers; }

  /* returns the view of all the superdrops in the domain */
  subviewd_supers domain_supers() const { return Kokkos::subview(totsupers, domainrefs); }

  /* returns the view of all the superdrops in the domain. read-only means superdrops in
  the view are const */
  subviewd_constsupers domain_supers_readonly() const {
    return Kokkos::subview(totsupers, domainrefs);
  }

  /* returns the total number of all the superdrops in the domain (excluding out of bounds ones) */
  size_t domain_nsupers() const { return domainrefs.second - domainrefs.first; }

  /* returns true if superdrops in view are sorted by their sdgbxindexes in ascending order */
  bool is_sorted() const { return sort_by_sdgbxindex.is_sorted(totsupers); }

  /* sort superdroplets by sdgbxindex and then (re-)set the totsupers view and the refs for the
  superdroplets that are within the domain (sdgbxindex within gbxindex_range for a given node) */
  viewd_supers sort_totsupers() {
    auto totsupers_ = sort_by_sdgbxindex(totsupers);
    set_totsupers_domainrefs(totsupers_);
    return totsupers;
  }

  /* Only use if you know what you're doing(!) Return leaves SupersInDomain instance in an
  intermediate state. Function sorts superdroplets by sdgbxindex but does not set the supers view
  nor the refs for the superdroplets that are within the domain. This means 'totsupers' may change,
  returned view may no longer be 'totsupers' and the domainrefs may be invalid. */
  viewd_supers sort_totsupers_without_set() { return sort_by_sdgbxindex(totsupers); }

  /* Only use if you know what you're doing(!) Assigns totsupers to given view and then sorts
  superdroplets by sdgbxindex with possible (re-)setting of the totsupers view and the refs for the
  superdroplets that are within the domain (sdgbxindex within gbxindex_range for a given node) */
  viewd_supers sort_and_set_totsupers(const viewd_supers totsupers_) {
    totsupers = totsupers_;
    return sort_totsupers();
  }
};

#endif  // LIBS_GRIDBOXES_SUPERSINDOMAIN_HPP_
