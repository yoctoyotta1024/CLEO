/*
 * ----- CLEO -----
 * File: findrefs.hpp
 * Project: gridboxes
 * Created Date: Wednesday 29th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 27th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions for findind references to superdrops
 * with particular sdgbxindex in a superdroplet view
 * (see it's use e.g. in supersingbx.hpp)
 */

#ifndef FINDREFS_HPP
#define FINDREFS_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "superdrops/kokkosaliases_sd.hpp"

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

template <typename ViewSupers>
inline kkpair_size_t find_refs(const ViewSupers totsupers,
                               unsigned int idx)
/* returns position in view of {first, last} superdrop
that occupies gridbox, ie. that has sdgbxindex == idx.
Function is outermost level of parallelism. */
{
  namespace SRP = SetRefPreds;
  const auto ref0 = size_t{find_ref(totsupers, SRP::Ref0{idx})};
  const auto ref1 = size_t{find_ref(totsupers, SRP::Ref1{idx})};

  return {ref0, ref1};
}

template <typename TeamMemberType, typename ViewSupers>
KOKKOS_INLINE_FUNCTION kkpair_size_t
find_refs(const TeamMemberType &team_member,
          const ViewSupers totsupers,
          unsigned int idx)
/* returns position in view of {first, last} superdrop
that occupies gridbox, ie. that has sdgbxindex == idx.
Function works within 1st layer of heirarchal
parallelism for a team_member of a league */
{
  namespace SRP = SetRefPreds;
  const auto ref0 = size_t{find_ref(team_member, totsupers, SRP::Ref0{idx})};
  const auto ref1 = size_t{find_ref(team_member, totsupers, SRP::Ref1{idx})};

  return {ref0, ref1};
}

template <typename Pred, typename ViewSupers>
inline size_t find_ref(const ViewSupers totsupers,
                       const Pred pred)
/* returns distance from begining of totsupers
view to the superdroplet that is first to fail
to satisfy given Predicate "pred".
Function is outermost level of parallelism. */
{
  namespace KE = Kokkos::Experimental;

  /* iterator to first superdrop in
  totsupers that fails to satisfy pred */
  const auto iter = KE::partition_point("find_ref",
                                        ExecSpace(),
                                        totsupers,
                                        pred);

  return makeref(totsupers, iter);
}

template <typename Pred, typename TeamMemberType, typename ViewSupers>
KOKKOS_INLINE_FUNCTION size_t
find_ref(const TeamMemberType &team_member,
         const ViewSupers totsupers,
         const Pred pred)
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
  const auto iter = KE::partition_point(team_member,
                                        totsupers,
                                        pred);
  return makeref(totsupers, iter);
}

template <typename Iter, typename ViewSupers>
KOKKOS_INLINE_FUNCTION size_t
makeref(const ViewSupers totsupers, const Iter iter)
/* makes ref (to use in refs pair for supers subview) by
returning distance from start of totsupers view to
position given by iterator 'iter'.
Note casting away signd-ness of distance. */
{
  namespace KE = Kokkos::Experimental;

  const auto ref = KE::distance(KE::begin(totsupers), iter);
  return static_cast<size_t>(ref);
}

template <typename ViewSupers>
inline kkpair_size_t find_domainrefs(const ViewSupers totsupers)
/* returns position in view of {first, last} superdrop
that is in domain, ie. that has sdgbxindex < oob_idx.
Function is outermost level of parallelism. */
{
  constexpr unsigned int oob_idx = LIMITVALUES::uintmax; // out of bounds sdgbxindex

  namespace SRP = SetRefPreds;
  const auto ref0 = size_t{0};
  const auto ref1 = size_t{find_ref(totsupers, SRP::Ref0{oob_idx})};

  return {ref0, ref1};
}

#endif // FINDREFS_HPP