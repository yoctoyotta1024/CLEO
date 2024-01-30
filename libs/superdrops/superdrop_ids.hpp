/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: superdrop_ids.hpp
 * Project: superdrops
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
 * File Description:
 * Header file for structs and functions for assignig
 * superdroplets with identifiers (IDs). E.g. and ID may
 * be a unique number starting from 0,
 * or non-existant (occupying no memory)
 */

#ifndef LIBS_SUPERDROPS_SUPERDROP_IDS_HPP_
#define LIBS_SUPERDROPS_SUPERDROP_IDS_HPP_

#include <ostream>

#include <Kokkos_Core.hpp>

/* struct containing value of SD identity (8bytes integer) */
struct IntID {
  size_t value;

  class Gen {
   public:
    /* note this generator is not thread safe
    (_idx++ undefined in multi-threaded environment) */
    IntID next() { return {_id++}; }

    /* note this generator assumes idx was
    thread safe generated (ie. is unique) */
    KOKKOS_INLINE_FUNCTION IntID next(const size_t id) { return {id}; }

   private:
    size_t _id = 0;
  };
};

/* struct for non-existant SD identity */
struct EmptyID {
  class Gen {
   public:
    KOKKOS_INLINE_FUNCTION EmptyID next() { return {}; }

    KOKKOS_INLINE_FUNCTION EmptyID next(const unsigned int kk) { return {}; }
  };
};

/* print SD identity given it is an IntID (an integer) */
inline std::ostream &operator<<(std::ostream &os, const IntID &id) {
  os << id.value;
  return os;
}

/* print null statement given SD identity is EmptyID (non-existent) */
inline std::ostream &operator<<(std::ostream &os, const EmptyID &id) {
  os << "(Undefined) No ID";
  return os;
}

#endif  // LIBS_SUPERDROPS_SUPERDROP_IDS_HPP_
