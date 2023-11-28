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
 * functionailty to handle superdroplets
 * inside CLEO's gridboxes
 */

#include "./supersingbx.hpp"

struct Pred
/* predicate to check superdrop
has matching sdgbxindex to ii*/
{
  unsigned int ii;

  KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
  {
    return op.get_sdgbxindex() == ii;
  }
};

bool SupersInGbx::iscorrect() const
/* assumes supers is already sorted via sdgbxindex. checks that all
superdrops in view which have matching sdgbxindex to gbxindex are indeed
included in (*this) subview (according to refs). Three criteria must
be true for iscorrect to return true: (1) all superdrops in current
subview have matching index. (2) all superdrops preceeding current
subview do not have matching index. (3) all superdrops after current
subview also do not have matching index. */
{
  const Pred pred{idx};
  const auto crit1(is_pred(pred));
  const auto crit2(is_prednot(pred, {0, refs.first}));
  const auto crit3(is_prednot(pred, {refs.second, totsupers.extent(0)}));

  return (crit1 && crit2 && crit3);
}

template <typename Pred>
bool SupersInGbx::is_pred(const Pred pred) const
/* returns true if all superdrops in subview
between refs satisfy the Predicate "pred" */
{
  return Kokkos::Experimental::
      all_of("is_pred",
             Kokkos::DefaultExecutionSpace(), // should match kokkosaliases.hpp
             (*this)(), pred);
}

template bool SupersInGbx::
    is_pred<Pred>(const Pred) const;

template <typename Pred>
bool SupersInGbx::is_prednot(const Pred pred,
                             const SupersInGbx::kkpair refs4pred) const
/* returns true if all superdrops in subview
between r0 and r1 do not satisfy pred */
{
  const subviewd_constsupers supers4pred(Kokkos::subview(totsupers,
                                                         refs4pred));

  return Kokkos::Experimental::
      none_of("is_prednot",
              Kokkos::DefaultExecutionSpace(), // should match kokkosaliases.hpp
              supers4pred, pred);
}

template bool SupersInGbx::
    is_prednot<Pred>(const Pred, const SupersInGbx::kkpair) const;

template <typename Pred>
size_t SupersInGbx::find_ref(const Pred pred) const
/* returns distance from begining of totsupers
view to the superdroplet that is first to fail
to satisfy given Predicate "pred" */
{
  namespace KE = Kokkos::Experimental;

  /* iterator to first superdrop in
  totsupers that fails to satisfy pred */
  const auto iter(KE::partition_point("findref",
                                      Kokkos::DefaultExecutionSpace(), // should match kokkosaliases.hpp
                                      totsupers, pred));

  /* distance form start of totsupers
  (casting away signd-ness)*/
  const auto ref0 = KE::distance(KE::begin(totsupers), iter);
  return static_cast<size_t>(ref0);
}

template size_t SupersInGbx::
    find_ref<SetRefPreds::Ref0>(const SetRefPreds::Ref0) const;

template size_t SupersInGbx::
    find_ref<SetRefPreds::Ref1>(const SetRefPreds::Ref1) const;