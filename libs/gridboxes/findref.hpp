/*
 * ----- CLEO -----
 * File: findref.hpp
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
 * functions for findind references to superdrops
 * with particular sdgbxindex in a superdroplet view
 * (see it's use e.g. in supersingbx.hpp)
 */

#ifndef FINDREF_HPP
#define FINDREF_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

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
  const auto iter(KE::partition_point("findref",
                                      ExecSpace(),
                                      totsupers,
                                      pred));

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
  const auto iter(KE::partition_point(team_member,
                                      totsupers,
                                      pred));
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

  const auto ref0 = KE::distance(KE::begin(totsupers), iter);
  return static_cast<size_t>(ref0);
}

#endif // FINDREF_HPP