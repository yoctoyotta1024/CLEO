/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: gbxindex.hpp
 * Project: gridboxes
 * Created Date: Tuesday 28th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions and structures related to
 * the unique indexes that label CLEO's gridboxes
 */

#ifndef LIBS_GRIDBOXES_GBXINDEX_HPP_
#define LIBS_GRIDBOXES_GBXINDEX_HPP_

/* struct containing gridbox index and its generator struct */
struct Gbxindex {
  unsigned int value;

  class Gen {
   public:
    /* note this generator is not thread safe
    (_idx++ undefined in multi-threaded environment) */
    Gbxindex next() { return {_idx++}; }

    /* note this generator assumes idx was
    thread safe generated (ie. is unique) */
    KOKKOS_INLINE_FUNCTION
    Gbxindex next(const unsigned int idx) { return {idx}; }

   private:
    unsigned int _idx = 0;
  };
};

#endif  // LIBS_GRIDBOXES_GBXINDEX_HPP_
