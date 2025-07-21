/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: supersingbx.hpp
 * Project: gridboxes
 * Created Date: Tuesday 28th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
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
KOKKOS_FUNCTION
bool SupersInGbx::iscorrect(const TeamMember &team_member,
                            const viewd_constsupers totsupers) const {
  const auto pred = Pred{idx};
  const auto crit1(is_pred(team_member, pred, totsupers));
  const auto crit2(is_prednot(team_member, pred, totsupers, {0, refs.first}));
  const auto crit3(is_prednot(team_member, pred, totsupers, {refs.second, totsupers.extent(0)}));

  return (crit1 && crit2 && crit3);
}

/* returns true if all superdrops in subview
between refs satisfy the Predicate "pred" */
template <typename Pred>
KOKKOS_FUNCTION bool SupersInGbx::is_pred(const TeamMember &team_member, const Pred pred,
                                          const viewd_constsupers totsupers) const {
  return Kokkos::Experimental::all_of(team_member, readonly(totsupers), pred);
}

/* returns true if all superdrops in subview between r0 and r1 do not satisfy pred */
template <typename Pred>
KOKKOS_FUNCTION bool SupersInGbx::is_prednot(const TeamMember &team_member, const Pred pred,
                                             const viewd_constsupers totsupers,
                                             const kkpair_size_t refs4pred) const {
  const auto supers4pred = Kokkos::subview(totsupers, refs4pred);
  return Kokkos::Experimental::none_of(team_member, supers4pred, pred);
}

template KOKKOS_FUNCTION bool SupersInGbx::is_pred<Pred>(const TeamMember &, const Pred,
                                                         const viewd_constsupers) const;

template KOKKOS_FUNCTION bool SupersInGbx::is_prednot<Pred>(const TeamMember &, const Pred,
                                                            const viewd_constsupers,
                                                            const kkpair_size_t) const;
