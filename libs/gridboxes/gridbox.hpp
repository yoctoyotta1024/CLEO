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
    supers_constview_type supers;               // reference to all superdrops view 
    unsigned int ii;                            // value of gbxindex which sdgbxindex of superdrops must match
    Kokkos::pair<size_t, size_t> refs = {0, 0}; // position in view of (first, last) superdrop that occupies gridbox
    
    using supers_subview = Kokkos::
        Subview<supers_constview_type, Kokkos::pair<size_t, size_t>>;

    template <typename Pred>
    KOKKOS_INLINE_FUNCTION size_t find_ref(Pred pred)
    /* returns distance from begining of supers view to
    the superdroplet that is first in supers to fail
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

  public:
    KOKKOS_INLINE_FUNCTION SupersInGbx() = default;  // Kokkos requirement for a (dual)View
    KOKKOS_INLINE_FUNCTION ~SupersInGbx() = default; // Kokkos requirement for a (dual)View

    KOKKOS_INLINE_FUNCTION SupersInGbx(supers_constview_type isupers,
                                       const unsigned int ii)
        : supers(isupers), ii(ii), refs(set_refs()) {}

    KOKKOS_INLINE_FUNCTION
    Kokkos::pair<size_t, size_t> set_refs()
    /* assumes supers is already sorted via sdgbxindex */
    {
      struct Ref0
      {
        unsigned int ii;

        KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
        {
          return op.get_sdgbxindex() < ii;
        }
      };

      struct Ref1
      {
        unsigned int ii;

        KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
        {
          return op.get_sdgbxindex() <= ii;
        }
      };

      return {find_ref(Ref0{ii}), find_ref(Ref1{ii})};
    }

    KOKKOS_INLINE_FUNCTION supers_subview operator()() const
    /* returns subview from view of superdrops
    refering to superdrops which occupy given gridbox */
    {
      return Kokkos::subview(supers, refs);
    }

    KOKKOS_INLINE_FUNCTION size_t nsupers() const
    {
      return refs.second - refs.first;
    }

    KOKKOS_INLINE_FUNCTION bool iscorrect() const
    /* assumes supers is already sorted via sdgbxindex */
    {
      std::cout << "\n--- set_refs() --- \n";
      for (size_t kk(0); kk < supers.extent(0); ++kk)
      {
        std::cout << supers(kk).id.value <<", ";
      }
      std::cout << "\n";
      for (size_t kk(0); kk < supers.extent(0); ++kk)
      {
        std::cout << supers(kk).get_sdgbxindex() <<", ";
      }
      std::cout << "--- --- --- --- ---\n";


      return true; //TODO: do this properly
      // see kokkos std algorithms all_of and none_of
    }
  };

public:
  struct Gbxindex
  /* struct containing gridbox index
  and its generator struct */
  {
    unsigned int value;
    class Gen
    {
    public:
      Gbxindex next() { return {_id++}; }

    private:
      unsigned int _id = 0;
    };
  };

  Gbxindex gbxindex;    // index (unique identifier) of gridbox
  State state;          // dynamical state of gridbox (e.g. thermodynamics)
  SupersInGbx supersingbx; // reference(s) to superdrops occupying gridbox
  Detectors detectors;  // detectors of various quantities

  KOKKOS_INLINE_FUNCTION Gridbox() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Gridbox() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  Gridbox(const Gbxindex igbxindex,
          const State istate,
          supers_constview_type supers)
      /* assumes supers view (or subview) already sorted via sdgbxindex */
      : gbxindex(igbxindex),
        state(istate),
        supersingbx(supers, gbxindex.value),
        detectors()
  {
  }

  KOKKOS_INLINE_FUNCTION
  auto get_gbxindex() {return gbxindex.value;}
};

#endif // GRIDBOX_HPP