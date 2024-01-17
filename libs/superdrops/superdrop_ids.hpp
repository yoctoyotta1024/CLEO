/*
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
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for structs and functions for assignig
 * superdroplets with identifiers (IDs). E.g. and ID may
 * be a unique number starting from 0,
 * or non-existant (occupying no memory)
 */

#ifndef SUPERDROP_IDS_HPP
#define SUPERDROP_IDS_HPP

#include <ostream>

#include <Kokkos_Core.hpp>

struct IntID
{
  /* struct containing value of SD identity (8bytes integer) */
  size_t value;

  class Gen
  {
  public:
    IntID next()
    /* note this generator is not thread safe
    (_idx++ undefined in multi-threaded environment) */
    {
      return {_id++};
    }

    KOKKOS_INLINE_FUNCTION IntID next(const size_t id)
    /* note this generator assumes idx was
    thread safe generated (ie. is unique) */
    {
      return {id};
    }

  private:
    size_t _id = 0;
  };
};

struct EmptyID
{
  /* struct for non-existant SD identity */
  class Gen
  {
  public:
    KOKKOS_INLINE_FUNCTION EmptyID next() { return {}; }

    KOKKOS_INLINE_FUNCTION EmptyID next(const unsigned int kk)
    {
      return {};
    }
  };
};

inline std::ostream &operator<<(std::ostream &os, const IntID &id)
/* print SD identity given it is an IntID (an integer) */
{
  os << id.value;
  return os;
}

inline std::ostream &operator<<(std::ostream &os, const EmptyID &id)
/* print null statement given SD identity is EmptyID (non-existent) */
{
  os << "(Undefined) No ID";
  return os;
}

#endif // SUPERDROP_IDS_HPP
