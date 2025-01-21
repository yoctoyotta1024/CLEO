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
#include <Kokkos_StdAlgorithms.hpp>

#include "gridboxes/findrefs.hpp"
#include "superdrops/kokkosaliases_sd.hpp"

/* References to identify the chunk of memory containing super-droplets occupying domain
(i.e. within any of the gridboxes on a single node), e.g. through std::span or Kokkos::subview) */
struct SupersInDomain {
 private:
  viewd_supers totsupers; /**< view of all superdrops (both in and out of bounds of domain) */

 public:
  SupersInDomain() = default;   // Kokkos requirement for a (dual)View
  ~SupersInDomain() = default;  // Kokkos requirement for a (dual)View

  explicit SupersInDomain(const viewd_supers i_totsupers) : totsupers(i_totsupers) {}

  void set_totsupers(const viewd_supers i_totsupers) { totsupers = i_totsupers; }

  viewd_supers get_totsupers() const { return totsupers; }

  /* returns the view of all the superdrops in the domain */
  subviewd_supers domain_supers() const {
    const auto domainrefs = find_domainrefs(totsupers);
    return Kokkos::subview(totsupers, domainrefs);
  }

  /* returns the view of all the superdrops in the domain. read-only means superdrops in
  the view are const */
  subviewd_constsupers domain_totsupers_readonly() const {
    const auto domainrefs = find_domainrefs(totsupers);
    return Kokkos::subview(totsupers, domainrefs);
  }

  /* returns the total number of all the superdrops in the domain (excluding out of bounds ones) */
  size_t domain_nsupers() const {
    const auto domainrefs = find_domainrefs(totsupers);
    return domainrefs.second - domainrefs.first;
  }
};

#endif  // LIBS_GRIDBOXES_SUPERSINDOMAIN_HPP_
