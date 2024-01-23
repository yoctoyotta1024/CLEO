/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
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
 * File Description:
 * functionailty to handle superdroplets
 * inside CLEO's gridboxes
 */

#include "gridboxes/supersingbx.hpp"

/* predicate to check superdrop
has matching sdgbxindex to ii*/
struct Pred {
  unsigned int ii;

  KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const {
    return op.get_sdgbxindex() == ii;
  }
};

/* assumes supers is already sorted via sdgbxindex. checks that all
superdrops in view which have matching sdgbxindex to gbxindex are indeed
included in (*this) subview (according to refs). Three criteria must
be true for iscorrect to return true: (1) all superdrops in current
subview have matching index. (2) all superdrops preceeding current
subview do not have matching index. (3) all superdrops after current
subview also do not have matching index. */
bool SupersInGbx::iscorrect() const {
  const Pred pred{idx};
  const auto crit1(is_pred(pred));
  const auto crit2(is_prednot(pred, {0, refs.first}));
  const auto crit3(is_prednot(pred, {refs.second, totsupers.extent(0)}));

  return (crit1 && crit2 && crit3);
}

/* returns true if all superdrops in subview
between refs satisfy the Predicate "pred" */
template <typename Pred>
bool SupersInGbx::is_pred(const Pred pred) const {
  namespace KE = Kokkos::Experimental;

  return KE::all_of("is_pred", ExecSpace(), (*this)(), pred);
}

/* returns true if all superdrops in subview
between r0 and r1 do not satisfy pred */
template <typename Pred>
bool SupersInGbx::is_prednot(const Pred pred, const kkpair_size_t refs4pred) const {
  namespace KE = Kokkos::Experimental;

  const subviewd_constsupers supers4pred(Kokkos::subview(totsupers, refs4pred));

  return KE::none_of("is_prednot", ExecSpace(), supers4pred, pred);
}

template bool SupersInGbx::is_pred<Pred>(const Pred) const;

template bool SupersInGbx::is_prednot<Pred>(const Pred, const kkpair_size_t) const;
