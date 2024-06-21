/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: supersingbx.hpp
 * Project: gridboxes
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions and structures related to
 * handling superdroplets inside CLEO's gridboxes
 */

#ifndef LIBS_GRIDBOXES_SUPERSINGBX_HPP_
#define LIBS_GRIDBOXES_SUPERSINGBX_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "gridboxes/findrefs.hpp"
#include "superdrops/kokkosaliases_sd.hpp"

/* References to identify the chunk of memory
containing super-droplets occupying a given Gridbox
(e.g. through std::span or Kokkos::subview) */
struct SupersInGbx {
 private:
  viewd_supers
      totsupers; /**< reference to view of all superdrops (both in and out of bounds of domain) */
  unsigned int idx;   /**< value of gbxindex which sdgbxindex of superdrops must match */
  kkpair_size_t refs; /**< position in view of (first, last) superdrop that occupies gridbox */

  template <typename Pred>
  bool is_pred(const Pred pred) const;
  /* returns true if all superdrops in subview
  between refs satisfy the Predicate "pred" */

  template <typename Pred>
  bool is_prednot(const Pred pred, const kkpair_size_t refs4pred) const;
  /* returns true if all superdrops in subview
  between r0 and r1 do not satisfy pred */

 public:
  SupersInGbx() = default;   // Kokkos requirement for a (dual)View
  ~SupersInGbx() = default;  // Kokkos requirement for a (dual)View

  /* assumes supers view (or subview) already sorted via sdgbxindex. Constructor
  works outside of parallelism */
  SupersInGbx(const viewd_supers i_totsupers, const unsigned int i_idx)
      : totsupers(i_totsupers), idx(i_idx), refs({0, 0}) {
    set_refs();
  }

  /* assumes supers view (or subview) already sorted via sdgbxindex. Constructor
  works within parallel team policy on host given member 'team_member' */
  SupersInGbx(const viewd_supers i_totsupers, const unsigned int i_idx, const kkpair_size_t i_refs)
      : totsupers(i_totsupers), idx(i_idx), refs(i_refs) {}

  /* assumes totsupers is already sorted via sdgbxindex. checks that all
  superdrops in view which have matching sdgbxindex to idx are indeed
  included in (*this) subview (according to refs). Three criteria must
  be true for iscorrect to return true: (1) all superdrops in current
  subview have matching index. (2) all superdrops preceeding current
  subview do not have matching index. (3) all superdrops after current
  subview also do not have matching index. */
  bool iscorrect() const;

  /* assumes totsupers is already sorted via sdgbxindex.
  sets 'refs' to pair with positions of first and last
  superdrops in view which have matching sdgbxindex to idx.
  Function is outside of parallelism (ie. in serial code). */
  inline void set_refs();

  /* assumes totsupers is already sorted via sdgbxindex.
  sets 'refs' to pair with positions of first and last
  superdrops in view which have matching sdgbxindex to idx.
  Function works within 1st layer of heirarchal parallelism
  for a team_member of a league */
  KOKKOS_INLINE_FUNCTION void set_refs(const TeamMember &team_member);

  /* returns subview from view of superdrops referencing superdrops
  which occupy given gridbox (according to refs) */
  KOKKOS_INLINE_FUNCTION
  subviewd_supers operator()() const { return Kokkos::subview(totsupers, refs); }

  /* returns read-only subview from view of superdrops
  referencing superdrops which occupy given gridbox
  (according to refs). read-only means superdrops
  in the subview are const */
  KOKKOS_INLINE_FUNCTION
  subviewd_constsupers readonly() const { return Kokkos::subview(totsupers, refs); }

  /* returns current number of superdrops referred to by gridbox */
  KOKKOS_INLINE_FUNCTION size_t nsupers() const { return refs.second - refs.first; }

  /* returns the view of all the superdrops in the domain.
  read-only means superdrops in the view are const */
  viewd_constsupers domain_totsupers_readonly() const {
    const auto domainrefs = find_domainrefs(totsupers);
    return Kokkos::subview(totsupers, domainrefs);
  }

  /* returns the total number of all the superdrops in the domain */
  size_t domain_totnsupers() const {
    const auto domainrefs = find_domainrefs(totsupers);
    return domainrefs.second - domainrefs.first;
  }
};

/* assumes totsupers is already sorted via sdgbxindex.
sets 'refs' to pair with positions of first and last
superdrops in view which have matching sdgbxindex to idx.
Function is outside of parallelism (ie. in serial code). */
inline void SupersInGbx::set_refs() { refs = find_refs(totsupers, idx); }

/* assumes totsupers is already sorted via sdgbxindex.
sets 'refs' to pair with positions of first and last
superdrops in view which have matching sdgbxindex to idx.
Function works within 1st layer of heirarchal
parallelism for a team_member of a league */
KOKKOS_INLINE_FUNCTION void SupersInGbx::set_refs(const TeamMember &team_member) {
  const auto new_refs = find_refs(team_member, totsupers, idx);

  Kokkos::single(
      Kokkos::PerTeam(team_member), [new_refs](kkpair_size_t &refs) { refs = new_refs; }, refs);
}

#endif  // LIBS_GRIDBOXES_SUPERSINGBX_HPP_
