/*
 * ----- CLEO -----
 * File: gridbox.cpp
 * Project: gridboxes
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 19th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Functions and structures related to the CLEO gridboxes
 */

#include "./gridbox.hpp"

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

bool Gridbox::SupersInGbx::iscorrect() const
/* assumes supers is already sorted via sdgbxindex. checks that all
superdrops in view which have matching sdgbxindex to ii are indeed
included in (*this) subview (according to refs). Three criteria must
be true for iscorrect to return true: (1) all superdrops in current
subview have matching index. (2) all superdrops preceeding current
subview do not have matching index. (3) all superdrops after current
subview also do not have matching index. */
{
  const Pred pred{ii};
  const auto crit1(is_pred(pred));
  const auto crit2(is_prednot(pred, {0, refs.first}));
  const auto crit3(is_prednot(pred, {refs.second, supers.extent(0)}));

  return (crit1 && crit2 && crit3);
}