// Author: Clara Bayley
// File: superdrop_ids.hpp
/* Header file for structs and functions for assignig
superdroplets with identifiers (IDs). E.g. and ID may 
be a unique number starting from 0, 
or non-existant (occupying no memory) */ 

#ifndef SUPERDROP_IDS_HPP
#define SUPERDROP_IDS_HPP

#include <iostream>

struct IntID {
  /* struct containing value of SD identity (8bytes integer) */
  size_t value;

  class Gen {
  public:
    IntID next () { return {_id++}; }
  private:
    size_t _id = 0;
  };
};

struct EmptyID {
  /* struct for non-existant SD identity */
  class Gen {
  public:
    EmptyID next () { return {}; }
  };
};

inline std::ostream& operator << (std::ostream& os, const IntID &id)
/* print SD identity given it is an IntID (an integer) */
{
   os << id.value;
   return os;
}

inline std::ostream& operator << (std::ostream& os, const EmptyID &id)
/* print null statement given SD identity is EmptyID (non-existent) */
{
  os << "(Undefined) No ID";
  return os;
}

#endif // SUPERDROP_IDS_HPP