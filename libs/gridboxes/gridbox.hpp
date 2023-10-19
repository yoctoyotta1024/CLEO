/*
 * ----- CLEO -----
 * File: gridbox.hpp
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

#ifndef GRIDBOX_HPP
#define GRIDBOX_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "./detectors.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

struct Gridbox
/* each gridbox has unique identifier and contains a
reference to superdroplets in gridbox, alongside the
Gridbox's State (e.g. thermodynamic variables
used for SDM) and detectors for tracking chosen variables */
{
private:
  using supers_constview_type = Kokkos::View<const Superdrop *>; // should match that in kokkosaliases.hpp
  
  struct SupersInGbx
  /* References to identify the chunk of memory
  containing super-droplets occupying a given Gridbox
  (e.g. through std::span or Kokkos::subview) */
  {
  private:
    using supers_subview = Kokkos::Subview<supers_constview_type,
                                           Kokkos::pair<size_t, size_t>>;
    supers_constview_type supers; // reference to all superdrops view
    unsigned int ii;              // value of gbxindex which sdgbxindex of superdrops must match
    Kokkos::pair<size_t, size_t> refs = {0, 0}; // position in view of (first, last) superdrop that occupies gridbox

    template <typename Pred>
    inline size_t find_ref(const Pred pred) const;
    /* returns distance from begining of supers view
    to the superdroplet that is first to fail
    to satisfy given Predicate "pred" */

    template <typename Pred>
    inline bool is_pred(const Pred pred) const;
    /* returns true if all superdrops in subview
    between refs satisfy the Predicate "pred" */

    template <typename Pred>
    inline bool is_prednot(const Pred pred,
                    const Kokkos::pair<size_t, size_t> refs4pred) const;
    /* returns true if all superdrops in subview
    between r0 and r1 do not satisfy pred */

  public:
    KOKKOS_INLINE_FUNCTION SupersInGbx() = default;  // Kokkos requirement for a (dual)View
    KOKKOS_INLINE_FUNCTION ~SupersInGbx() = default; // Kokkos requirement for a (dual)View

    SupersInGbx(const supers_constview_type isupers,
                const unsigned int ii)
        : supers(isupers), ii(ii), refs(set_refs()) {}

    inline Kokkos::pair<size_t, size_t> set_refs();
    /* assumes supers is already sorted via sdgbxindex.
    returns pair which are positions of first and last
    superdrops in view which have matching sdgbxindex to ii */

    bool iscorrect() const;
    /* assumes supers is already sorted via sdgbxindex. checks that all 
    superdrops in view which have matching sdgbxindex to ii are indeed
    included in (*this) subview (according to refs). Three criteria must
    be true for iscorrect to return true: (1) all superdrops in current
    subview have matching index. (2) all superdrops preceeding current
    subview do not have matching index. (3) all superdrops after current
    subview also do not have matching index. */

    KOKKOS_INLINE_FUNCTION
    supers_subview operator()() const
    /* returns subview from view of superdrops referencing superdrops
    which occupy given gridbox (according to refs) */
    {
      return Kokkos::subview(supers, refs);
    }

    KOKKOS_INLINE_FUNCTION size_t nsupers() const
    /* returns current number of superdrops referred to by gridbox */
    {
      return refs.second - refs.first;
    }
  };

public:
  struct Gbxindex
  /* struct containing gridbox index and its generator struct */
  {
    unsigned int value;
    class Gen
    {
    public:
      KOKKOS_INLINE_FUNCTION Gbxindex next() { return {_id++}; }

    private:
      unsigned int _id = 0;
    };
  };

  Gbxindex gbxindex;       // index (unique identifier) of gridbox
  State state;             // dynamical state of gridbox (e.g. thermodynamics)
  SupersInGbx supersingbx; // reference(s) to superdrops occupying gridbox
  Detectors detectors;     // detectors of various quantities

  KOKKOS_INLINE_FUNCTION Gridbox() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Gridbox() = default; // Kokkos requirement for a (dual)View

  Gridbox(const Gbxindex igbxindex,
          const State istate,
          const supers_constview_type supers)
      /* assumes supers view (or subview) already sorted via sdgbxindex */
      : gbxindex(igbxindex),
        state(istate),
        supersingbx(supers, gbxindex.value),
        detectors()
  {
  }

  KOKKOS_INLINE_FUNCTION
  auto get_gbxindex() const { return gbxindex.value; }
};

template <typename Pred>
inline size_t Gridbox::SupersInGbx::find_ref(const Pred pred) const
/* returns distance from begining of supers view
to the superdroplet that is first to fail
to satisfy given Predicate "pred" */
{
  namespace KE = Kokkos::Experimental;

  /* iterator to first superdrop in
  supers that fails to satisfy pred */
  const auto iter(KE::partition_point("findref",
                                      Kokkos::DefaultExecutionSpace(),
                                      supers, pred));

  /* distance form start of supers
  (casting away signd-ness)*/
  const auto ref0 = KE::distance(KE::begin(supers), iter);
  return static_cast<size_t>(ref0);
}

template <typename Pred>
inline bool Gridbox::SupersInGbx::is_pred(const Pred pred) const
/* returns true if all superdrops in subview
between refs satisfy the Predicate "pred" */
{
  return Kokkos::Experimental::
      all_of("is_pred",
             Kokkos::DefaultExecutionSpace(),
             (*this)(), pred);
}

template <typename Pred>
inline bool Gridbox::SupersInGbx::
    is_prednot(const Pred pred,
               const Kokkos::pair<size_t, size_t> refs4pred) const
/* returns true if all superdrops in subview
between r0 and r1 do not satisfy pred */
{
  const supers_subview
      supers4pred(Kokkos::subview(supers, refs4pred));

  return Kokkos::Experimental::
      none_of("is_prednot",
              Kokkos::DefaultExecutionSpace(),
              supers4pred, pred);
}

inline Kokkos::pair<size_t, size_t>
Gridbox::SupersInGbx::set_refs()
/* assumes supers is already sorted via sdgbxindex.
returns pair which are positions of first and last
superdrops in view which have matching sdgbxindex to ii */
{
  struct Ref0
  /* predicate to find first superdrop in
  view which has matching sdgbxindex to ii */
  {
    unsigned int ii;

    KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
    {
      return op.get_sdgbxindex() < ii;
    }
  };

  struct Ref1
  /* predicate to find last superdrop in
  view which has matching sdgbxindex to ii */
  {
    unsigned int ii;

    KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
    {
      return op.get_sdgbxindex() <= ii;
    }
  };

  return {find_ref(Ref0{ii}), find_ref(Ref1{ii})};
}

#endif // GRIDBOX_HPP