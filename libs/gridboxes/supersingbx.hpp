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

/* References to identify the chunk of memory containing super-droplets occupying a given Gridbox.
You must self ensure std::span/Kokkos::(sub)View used to find super-droplets is correct for
current refs */
struct SupersInGbx {
 private:
  unsigned int idx;   /**< value of gbxindex which sdgbxindex of superdrops must match */
  kkpair_size_t refs; /**< position in view of (first, last) superdrop that occupies gridbox */

  /* returns true if all superdrops in subview between refs satisfy the Predicate "pred" */
  template <typename Pred>
  KOKKOS_FUNCTION bool is_pred(const TeamMember &team_member, const Pred pred,
                               const viewd_constsupers totsupers) const;

  /* returns true if all superdrops in subview between r0 and r1 do not satisfy pred */
  template <typename Pred>
  KOKKOS_FUNCTION bool is_prednot(const TeamMember &team_member, const Pred pred,
                                  const viewd_constsupers totsupers,
                                  const kkpair_size_t refs4pred) const;

 public:
  SupersInGbx() = default;   // Kokkos requirement for a (dual)View
  ~SupersInGbx() = default;  // Kokkos requirement for a (dual)View

  /* assumes supers view (or subview) already sorted via sdgbxindex. Constructor
  works outside of parallelism to find refs given sorted superdrops in domain */
  SupersInGbx(const unsigned int idx_, const subviewd_constsupers domainsupers)
      : idx(idx_), refs({0, 0}) {
    set_refs(domainsupers);
  }

  /* Constructor works within parallel team policy on host given member 'team_member' */
  SupersInGbx(const unsigned int idx_, const kkpair_size_t refs_) : idx(idx_), refs(refs_) {}

  /* assumes domainsupers is already sorted via sdgbxindex. checks that all
  superdrops in view which have matching sdgbxindex to idx are indeed
  included in (*this) subview (according to refs). Three criteria must
  be true for iscorrect to return true: (1) all superdrops in current
  subview have matching index. (2) all superdrops preceeding current
  subview do not have matching index. (3) all superdrops after current
  subview also do not have matching index. */
  KOKKOS_FUNCTION bool iscorrect(const TeamMember &team_member,
                                 const viewd_constsupers totsupers) const;

  /* assumes domainsupers is already sorted via sdgbxindex.
  sets 'refs' to pair with positions of first and last
  superdrops in view which have matching sdgbxindex to idx.
  Function is outside of parallelism (ie. in serial code). */
  KOKKOS_INLINE_FUNCTION void set_refs(const subviewd_constsupers domainsupers);

  /* assumes domainsupers is already sorted via sdgbxindex.
  sets 'refs' to pair with positions of first and last
  superdrops in view which have matching sdgbxindex to idx.
  Function works within 1st layer of heirarchal parallelism
  for a team_member of a league */
  KOKKOS_INLINE_FUNCTION void set_refs(const TeamMember &team_member,
                                       const subviewd_constsupers domainsupers);

  /* returns subview from view of superdrops referencing superdrops
  which occupy given gridbox (according to refs) */
  KOKKOS_INLINE_FUNCTION
  subviewd_supers operator()(const subviewd_supers domainsupers) const {
    return Kokkos::subview(domainsupers, refs);
  }

  /* returns read-only subview from view of superdrops referencing superdrops which
  occupy given gridbox (according to refs). read-only means superdrops in the subview are const */
  KOKKOS_INLINE_FUNCTION
  subviewd_constsupers readonly(const subviewd_constsupers domainsupers) const {
    return Kokkos::subview(domainsupers, refs);
  }

  /* returns current number of superdrops referred to by gridbox */
  KOKKOS_INLINE_FUNCTION size_t nsupers() const { return refs.second - refs.first; }
};

/* assumes domainsupers is already sorted via sdgbxindex.
sets 'refs' to pair with positions of first and last
superdrops in view which have matching sdgbxindex to idx.
Function is outside of parallelism (ie. in serial code). */
KOKKOS_INLINE_FUNCTION void SupersInGbx::set_refs(const subviewd_constsupers domainsupers) {
  refs = find_refs(domainsupers, idx);
}

/* assumes domainsupers is already sorted via sdgbxindex.
sets 'refs' to pair with positions of first and last
superdrops in view which have matching sdgbxindex to idx.
Function works within 1st layer of heirarchal
parallelism for a team_member of a league */
KOKKOS_INLINE_FUNCTION void SupersInGbx::set_refs(const TeamMember &team_member,
                                                  const subviewd_constsupers domainsupers) {
  const auto new_refs = find_refs(team_member, domainsupers, idx);
  Kokkos::single(
      Kokkos::PerTeam(team_member), [new_refs](kkpair_size_t &refs) { refs = new_refs; }, refs);
}

#endif  // LIBS_GRIDBOXES_SUPERSINGBX_HPP_
