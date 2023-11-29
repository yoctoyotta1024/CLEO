/*
 * ----- CLEO -----
 * File: supersingbx.hpp
 * Project: gridboxes
 * Created Date: Tuesday 28th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 28th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Functions and structures related to
 * handling superdroplets inside CLEO's gridboxes
 */

#ifndef SUPERSINGBX_HPP
#define SUPERSINGBX_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "superdrops/kokkosaliases_sd.hpp"

struct SupersInGbx
/* References to identify the chunk of memory
containing super-droplets occupying a given Gridbox
(e.g. through std::span or Kokkos::subview) */
{
private:
  using kkpair = Kokkos::pair<size_t, size_t>;

  viewd_supers totsupers; // reference to view of all superdrops (in total domain)
  unsigned int idx;       // value of gbxindex which sdgbxindex of superdrops must match
  kkpair refs;            // position in view of (first, last) superdrop that occupies gridbox

  template <typename Pred>
  size_t find_ref(const Pred pred) const;
  /* returns distance from begining of totsupers
  view to the superdroplet that is first to fail
  to satisfy given Predicate "pred".
  Function is outermost level of parallelism. */

  template <typename Pred, typename TeamMemberType>
  KOKKOS_INLINE_FUNCTION size_t
  find_ref(const TeamMemberType &team_member,
           const Pred pred) const;
  /* returns distance from begining of totsupers
  view to the superdroplet that is first to fail
  to satisfy given Predicate "pred".
  Function is 2nd level of nested parallelism,
  i.e. is thread parallelism within a league 
  for a given team_member */

  template <typename Pred>
  bool is_pred(const Pred pred) const;
  /* returns true if all superdrops in subview
  between refs satisfy the Predicate "pred" */

  template <typename Pred>
  bool is_prednot(const Pred pred,
                  const kkpair refs4pred) const;
  /* returns true if all superdrops in subview
  between r0 and r1 do not satisfy pred */

  template <typename Iter>
  KOKKOS_INLINE_FUNCTION
  size_t makeref(const Iter iter) const
  /* makes ref (to use in refs pair for supers subview) by
  returning distance from start of totsupers view to
  position given by iterator 'iter'.
  Note casting away signd-ness of distance. */
  {
    namespace KE = Kokkos::Experimental;
    
    const auto ref0 = KE::distance(KE::begin(totsupers), iter);
    return static_cast<size_t>(ref0);
  }

public:
  KOKKOS_INLINE_FUNCTION SupersInGbx() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~SupersInGbx() = default; // Kokkos requirement for a (dual)View

  SupersInGbx(const viewd_supers i_totsupers,
              const unsigned int i_idx)
      : totsupers(i_totsupers), idx(i_idx), refs({0, 0})
  /* assumes supers view (or subview) already
  sorted via sdgbxindex. Constructor works
  outside of parallelism */
  {
    set_refs();
  }

  SupersInGbx(const HostTeamMember &team_member,
              const viewd_supers i_totsupers,
              const unsigned int i_idx)
      : totsupers(i_totsupers), idx(i_idx), refs({0, 0})
  /* assumes supers view (or subview) already sorted
  via sdgbxindex. Constructor works within parallel
  team policy on host given member 'team_member' */
  {
    set_refs(team_member);
  }

  bool iscorrect() const;
  /* assumes totsupers is already sorted via sdgbxindex. checks that all
  superdrops in view which have matching sdgbxindex to idx are indeed
  included in (*this) subview (according to refs). Three criteria must
  be true for iscorrect to return true: (1) all superdrops in current
  subview have matching index. (2) all superdrops preceeding current
  subview do not have matching index. (3) all superdrops after current
  subview also do not have matching index. */

  inline void set_refs();
  /* assumes totsupers is already sorted via sdgbxindex.
  sets 'refs' to pair with positions of first and last
  superdrops in view which have matching sdgbxindex to idx.
  Function is outside of parallelism (ie. in serial code). */

  template <typename TeamMemberType>
  KOKKOS_INLINE_FUNCTION void
  set_refs(const TeamMemberType &team_member);
  /* assumes totsupers is already sorted via sdgbxindex.
  sets 'refs' to pair with positions of first and last
  superdrops in view which have matching sdgbxindex to idx.
  Function works within 1st layer of heirarchal parallelism
  for a team_member of a league */

  KOKKOS_INLINE_FUNCTION
  subviewd_supers operator()() const
  /* returns subview from view of superdrops referencing superdrops
  which occupy given gridbox (according to refs) */
  {
    return Kokkos::subview(totsupers, refs);
  }

  KOKKOS_INLINE_FUNCTION
  subviewd_constsupers readonly() const
  /* returns subview from view of superdrops referencing superdrops
  which occupy given gridbox (according to refs) */
  {
    return Kokkos::subview(totsupers, refs);
  }

  mirrorh_constsupers hostcopy() const
  /* returns mirror view on host for const supers in
  gridbox. If supers view is on device memory, a
  deep copy is performed */
  {
    const subviewd_constsupers d_supers = readonly();
    auto h_supers = Kokkos::create_mirror_view(d_supers);
    Kokkos::deep_copy(h_supers, d_supers);

    return h_supers;
  }

  KOKKOS_INLINE_FUNCTION size_t nsupers() const
  /* returns current number of superdrops referred to by gridbox */
  {
    return refs.second - refs.first;
  }

  KOKKOS_INLINE_FUNCTION size_t domain_totnsupers() const
  /* returns current total number of superdrops in domain */
  {
    return totsupers.extent(0);
  }

  KOKKOS_INLINE_FUNCTION viewd_supers domain_totsupers() const
  /* returns current total number of superdrops in domain */
  {
    return totsupers;
  }
};

namespace SetRefPreds
/* namespace containing values of
constants with dimensions */
{

  struct Ref0
  /* struct for SupersInGbx::set_refs()
  predicate to find first superdrop in
  view which has matching sdgbxindex to idx */
  {
    unsigned int idx;

    KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
    {
      return op.get_sdgbxindex() < idx;
    }
  };

  struct Ref1
  /* struct for SupersInGbx::set_refs()
  predicate to find last superdrop in
  view which has matching sdgbxindex to idx */
  {
    unsigned int idx;

    KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
    {
      return op.get_sdgbxindex() <= idx;
    }
  };
}

inline void SupersInGbx::set_refs()
/* assumes totsupers is already sorted via sdgbxindex.
sets 'refs' to pair with positions of first and last
superdrops in view which have matching sdgbxindex to idx.
Function is outside of parallelism (ie. in serial code). */
{
  namespace SRP = SetRefPreds;
  refs = {find_ref(SRP::Ref0{idx}),
          find_ref(SRP::Ref1{idx})};
}

template <typename TeamMemberType>
KOKKOS_INLINE_FUNCTION void
SupersInGbx::set_refs(const TeamMemberType &team_member)
/* assumes totsupers is already sorted via sdgbxindex.
sets 'refs' to pair with positions of first and last
superdrops in view which have matching sdgbxindex to idx.
Function works within 1st layer of heirarchal parallelism
for a team_member of a league */
{
  namespace SRP = SetRefPreds;
  const kkpair new_refs = {find_ref(team_member, SRP::Ref0{idx}),
                           find_ref(team_member, SRP::Ref1{idx})};

  Kokkos::single(
      Kokkos::PerTeam(team_member),
      [new_refs](kkpair &refs)
      { refs = new_refs; },
      refs);
}

template <typename Pred, typename TeamMemberType>
KOKKOS_INLINE_FUNCTION size_t
SupersInGbx::find_ref(const TeamMemberType &team_member,
                      const Pred pred) const
/* returns distance from begining of totsupers
view to the superdroplet that is first to fail
to satisfy given Predicate "pred".
Function is 2nd level of nested parallelism,
i.e. is thread parallelism within a league
for a given team_member */
{
  namespace KE = Kokkos::Experimental;

  /* iterator to first superdrop in
  totsupers that fails to satisfy pred */
  const auto iter(KE::partition_point(team_member,
                                      totsupers,
                                      pred));
  return makeref(iter);                                   
}

#endif // SUPERSINGBX_HPP