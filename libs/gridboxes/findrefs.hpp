/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: findrefs.hpp
 * Project: gridboxes
 * Created Date: Wednesday 29th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions for findind references to superdrops
 * with particular sdgbxindex in a superdroplet view
 * (see it's use e.g. in supersingbx.hpp)
 */

#ifndef LIBS_GRIDBOXES_FINDREFS_HPP_
#define LIBS_GRIDBOXES_FINDREFS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "../cleoconstants.hpp"
#include "superdrops/kokkosaliases_sd.hpp"

/* namespace containing values of constants with dimensions */
namespace SetRefPreds {

/* struct for SupersInGbx::set_refs(...) predicate to find _first_ superdrop in
view which has matching sdgbxindex to idx */
struct Ref0 {
  unsigned int idx;

  KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const {
    return op.get_sdgbxindex() < idx;
  }
};

/* struct for SupersInGbx::set_refs(...) predicate to find _last_ superdrop in
view which has matching sdgbxindex to idx */
struct Ref1 {
  unsigned int idx;

  KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const {
    return op.get_sdgbxindex() <= idx;
  }
};
}  // namespace SetRefPreds

/* returns position in view of {first, last} superdrop that occupies gridbox,
ie. that has sdgbxindex == idx. Function is outermost level of parallelism. */
template <typename ViewSupers>
KOKKOS_INLINE_FUNCTION kkpair_size_t find_refs(const ViewSupers totsupers, unsigned int idx) {
  namespace SRP = SetRefPreds;
  const auto ref0 = size_t{find_ref(totsupers, SRP::Ref0{idx})};
  const auto ref1 = size_t{find_ref(totsupers, SRP::Ref1{idx})};

  return {ref0, ref1};
}

/* returns position in view of {first, last} superdrop that occupies gridbox,
ie. that has sdgbxindex == idx. Function works within 1st layer of heirarchal
parallelism for a team_member of a league */
template <typename TeamMemberType, typename ViewSupers>
KOKKOS_INLINE_FUNCTION kkpair_size_t find_refs(const TeamMemberType &team_member,
                                               const ViewSupers totsupers, unsigned int idx) {
  namespace SRP = SetRefPreds;
  const auto ref0 = size_t{find_ref(team_member, totsupers, SRP::Ref0{idx})};
  const auto ref1 = size_t{find_ref(team_member, totsupers, SRP::Ref1{idx})};

  return {ref0, ref1};
}

/* returns distance from begining of totsupers view to the superdroplet that is first to fail
to satisfy given Predicate "pred". Function can be outside or inside first-level of parallelism.
Parallel equivalent outside parallelism (KOKKOS_INLINE_FUNCTION ->inline)
is experimental and appears to be slower (!):
```
namespace KE = Kokkos::Experimental;
const auto iter = KE::partition_point("find_ref", ExecSpace(), totsupers, pred);
return makeref(KE::begin(totsupers), iter);
```
*/
template <typename Pred, typename ViewSupers>
KOKKOS_INLINE_FUNCTION size_t find_ref(const ViewSupers totsupers, const Pred pred) {
  return find_partition_point(totsupers, pred, 0, totsupers.extent(0));
}

/* returns element access index from begining of totsupers view to the superdroplet that
is first to fail to satisfy given Predicate "pred". Function is 2nd level of nested parallelism,
i.e. is thread parallelism within a league for a given team_member. Parallel equivalent
is experimental and appears to be slower (!):
```
namespace KE = Kokkos::Experimental;
const auto start = KE::begin(totsupers);
const auto end = KE::end(totsupers);
const auto iter = KE::partition_point(team_member, start, end, pred);
return makeref(start, iter);
```
*/
template <typename Pred, typename TeamMemberType, typename ViewSupers>
KOKKOS_INLINE_FUNCTION size_t find_ref(const TeamMemberType &team_member,
                                       const ViewSupers totsupers, const Pred pred) {
  return find_partition_point(totsupers, pred, 0, totsupers.extent(0));
}

/* returns element access index from begining of totsupers view to the superdroplet that is
first to fail to satisfy given Predicate "pred". Function is Kokkos GPU compatible version
of std::partition_point specific to kokkos view types */
template <typename ViewSupers, typename Pred>
KOKKOS_INLINE_FUNCTION size_t find_partition_point(const ViewSupers totsupers, const Pred pred,
                                                   size_t first, size_t start_length) {
  for (auto length = start_length; 0 < length;) {
    size_t half = length / 2;
    size_t middle = first + half;
    if (pred(totsupers(middle))) {
      first = middle + 1;
      length -= (half + 1);
    } else
      length = half;
  }

  return first;
}

/* makes ref (to use in refs pair for supers subview) by returning distance from
first iterator (e.g. start of totsupers view) to position given by iterator 'iter'.
Note casting away signed-ness of distance. */
template <typename Iter>
KOKKOS_INLINE_FUNCTION size_t makeref(const Iter start, const Iter iter) {
  namespace KE = Kokkos::Experimental;
  const auto ref = KE::distance(start, iter);
  return static_cast<size_t>(ref);
}

/* returns position in view of {first, last} superdrop that is in domain. First is
position of first superdrop with sdgbxindex >= gbxindex_range.first. Last is position
of last superdrop with sdgbxindex <= gbxindex_range.second.
Function valid in outermost level (outside) of parallelism on host. */
template <typename ViewSupers>
inline kkpair_size_t find_domainrefs(
    const ViewSupers totsupers, const Kokkos::pair<unsigned int, unsigned int> gbxindex_range) {
  namespace SRP = SetRefPreds;
  const auto ref0 = size_t{find_ref(totsupers, SRP::Ref0{gbxindex_range.first})};
  const auto ref1 = size_t{find_ref(totsupers, SRP::Ref1{gbxindex_range.second})};

  return {ref0, ref1};
}

/* returns position in view of {first, last} superdrop that is in domain, where first is assumed to
be at 0th position. Last is position of last superdrop with sdgbxindex < gbxindex_max.
Function valid in outermost level (outside) of parallelism default implementation valid on host. */
template <typename ViewSupers>
inline kkpair_size_t find_domainrefs(const Kokkos::HostSpace &ex, const ViewSupers totsupers,
                                     const unsigned int gbxindex_max) {
  namespace SRP = SetRefPreds;
  const auto ref1 = size_t{find_ref(totsupers, SRP::Ref1{gbxindex_max})};

  return {0, ref1};
}

/* returns position in view of {first, last} superdrop that is in domain, where first is assumed to
be at 0th position. Last is position of last superdrop with sdgbxindex < gbxindex_max.
Function valid in outermost level (outside) of parallelism valid for CUDA devices. */
template <class ExecutionSpace, typename ViewSupers>
inline kkpair_size_t find_domainrefs(const ExecutionSpace &ex, const ViewSupers totsupers,
                                     const unsigned int gbxindex_max) {
  namespace SRP = SetRefPreds;
  auto ref1 = Kokkos::View<size_t[1]>("domainref1");
  Kokkos::parallel_for(
      "find_domainrefs_cuda", 1, KOKKOS_LAMBDA(const int i) {
        ref1(0) = size_t{find_ref(totsupers, SRP::Ref1{gbxindex_max})};
      });
  auto h_ref1 = Kokkos::create_mirror_view_and_copy(HostSpace(), ref1);

  return {0, h_ref1(0)};
}

#endif  // LIBS_GRIDBOXES_FINDREFS_HPP_
