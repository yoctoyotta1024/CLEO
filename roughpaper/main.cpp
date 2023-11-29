/*
 * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 17th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * rough paper for checking small things
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <limits>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_NestedSort.hpp>


struct Superdrop
{
private:
  unsigned int sdgbxindex;
public:

KOKKOS_INLINE_FUNCTION
void set_sdgbxindex(const unsigned int ii){sdgbxindex = ii;}

KOKKOS_INLINE_FUNCTION
    unsigned int get_sdgbxindex() const
{
  return sdgbxindex;
}
};

namespace SetRefPreds
/* namespace containing values of
constants with dimensions */
{

  struct Ref0
  /* struct for SupersInGbx::set_refs()
  predicate to find first superdrop in
  view which has matching sdgbxindex to idx */
  {
    unsigned int idx;

    KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
    {
      return op.get_sdgbxindex() < idx;
    }
  };

  struct Ref1
  /* struct for SupersInGbx::set_refs()
  predicate to find last superdrop in
  view which has matching sdgbxindex to idx */
  {
    unsigned int idx;

    KOKKOS_INLINE_FUNCTION bool operator()(const Superdrop &op) const
    {
      return op.get_sdgbxindex() <= idx;
    }
  };
}

using viewd_supers = Kokkos::View<Superdrop *>;
using viewd_constsupers = Kokkos::View<const Superdrop *>; // view in device memory of const superdroplets
using ExecSpace = Kokkos::DefaultExecutionSpace;
using kkpair = Kokkos::pair<size_t, size_t>;

using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
using TeamMember = TeamPolicy::member_type;

template <typename Pred>
size_t find_ref_new(const TeamMember &team_member,
                    const Pred pred,
                    viewd_supers totsupers)
{
  namespace KE = Kokkos::Experimental;

  /* iterator to first superdrop in
  totsupers that fails to satisfy pred */
  const auto iter(KE::partition_point(team_member,
                  totsupers,
                  pred));

  /* distance form start of totsupers (casting away signed-ness) */
  const auto ref0 = KE::distance(KE::begin(totsupers), iter);
  return static_cast<size_t>(ref0);
}

kkpair set_refs_new(const TeamMember &team_member,
                    const unsigned int ii,
                    viewd_supers totsupers)
{
  const kkpair new_refs = {find_ref_new(team_member, SetRefPreds::Ref0{ii}, totsupers),
                           find_ref_new(team_member, SetRefPreds::Ref1{ii}, totsupers)};
  return new_refs;
}

void main_new(const size_t nsupers, const size_t ngbxs)
{
  std::cout << "\n ========== MAIN NEW ========== \n";
    viewd_supers supers("supers", nsupers);
    auto h_supers = Kokkos::create_mirror_view(supers);
    for (size_t kk(0); kk < nsupers; ++kk)
    {
      const unsigned int ii(kk/2);
      h_supers(kk) = Superdrop();
      h_supers(kk).set_sdgbxindex(ii);
      std::cout << "ii: " << h_supers(kk).get_sdgbxindex() << "\n";
    }
    Kokkos::deep_copy(supers, h_supers);

    Kokkos::View<kkpair *> viewd_refs("gbxs", ngbxs);
    Kokkos::parallel_for(
        "move_superdrops_in_domain",
        TeamPolicy(ngbxs, Kokkos::AUTO()),
        KOKKOS_LAMBDA(const TeamMember &team_member) {
          const int ii = team_member.league_rank();

          kkpair refs{0, 0};
          refs = set_refs_new(team_member, ii, supers);
          std::cout << "refs: " << refs.first << ", " << refs.second << "\n";

          viewd_refs(ii) = refs;
        });

    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      const auto subview = Kokkos::subview(supers, viewd_refs(ii));
      const size_t n(subview.extent(0));
      std::cout << "---- gbx: ii = " << ii << " -----\n"
                << "refs: " << viewd_refs(ii).first
                << ", " << viewd_refs(ii).second << "\n";

      for (size_t kk(0); kk < n; ++kk)
      {
        std::cout << "ii: " << subview(kk).get_sdgbxindex() << "\n"; 
      }
      std::cout << "n = " << n << "\n"; 
    }
  std::cout << "\n ========== ========== ========== \n";
}

void main_old(const size_t nsupers, const size_t ngbxs);

int main(int argc, char *argv[])
{
  const size_t nsupers(12);
  const size_t ngbxs(13);

  Kokkos::initialize(argc, argv);
  {
    main_old(nsupers, ngbxs);
    main_new(nsupers, ngbxs);
  }
  Kokkos::finalize();
}


/* --- old algorithm --- */
kkpair set_refs_old(const unsigned int ii,
                    viewd_constsupers totsupers);

template <typename Pred>
size_t find_ref_old(const Pred pred,
                    viewd_constsupers totsupers);

void main_old(const size_t nsupers, const size_t ngbxs)
{
  std::cout << "\n ========== MAIN OLD ========== \n";
    viewd_supers supers("supers", nsupers);
    auto h_supers = Kokkos::create_mirror_view(supers);
    for (size_t kk(0); kk < nsupers; ++kk)
    {
      const unsigned int ii(kk/3+1);
      h_supers(kk) = Superdrop();
      h_supers(kk).set_sdgbxindex(ii);
      std::cout << "ii: " << h_supers(kk).get_sdgbxindex() << "\n";
    }
    Kokkos::deep_copy(supers, h_supers);

    Kokkos::View<kkpair *> viewd_refs("gbxs", ngbxs); 
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      kkpair refs{0, 0};
      refs = set_refs_old(ii, supers);
      std::cout << "refs: " << refs.first << ", " << refs.second << "\n";
      
      viewd_refs(ii) = refs;
    }

    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      const auto subview = Kokkos::subview(supers, viewd_refs(ii));
      const size_t n(subview.extent(0));
      std::cout << "---- gbx: ii = " << ii << " -----\n"
                << "refs: " << viewd_refs(ii).first
                << ", " << viewd_refs(ii).second << "\n";

      for (size_t kk(0); kk < n; ++kk)
      {
        std::cout << "ii: " << subview(kk).get_sdgbxindex() << "\n"; 
      }
      std::cout << "n = " << n << "\n"; 
    }
  std::cout << "\n ========== ========== ========== \n";
}

kkpair set_refs_old(const unsigned int ii, viewd_constsupers totsupers)
{
  namespace SRP = SetRefPreds;
  const kkpair new_refs = {find_ref_old(SRP::Ref0{ii}, totsupers),
                           find_ref_old(SRP::Ref1{ii}, totsupers)};
  return new_refs;
}

template <typename Pred>
size_t find_ref_old(const Pred pred,
                    viewd_constsupers totsupers)
{
  namespace KE = Kokkos::Experimental;

  /* iterator to first superdrop in
  totsupers that fails to satisfy pred */
  const auto iter(KE::partition_point("findref",
                                      ExecSpace(),
                                      totsupers,
                                      pred));

  /* distance form start of totsupers
  (casting away signd-ness)*/
  const auto ref0 = KE::distance(KE::begin(totsupers), iter);
  return static_cast<size_t>(ref0);
}

template <class ExecutionSpace, class DataType, class... Properties>
void sort([[maybe_unused]] const ExecutionSpace& exec,
          const Kokkos::View<DataType, Properties...>& view) {...}

template <class DataType, class... Properties>
void sort(const Kokkos::View<DataType, Properties...>& view) {...}

template <class ExecutionSpace, class ComparatorType, class DataType,
          class... Properties>
void sort([[maybe_unused]] const ExecutionSpace& exec,
          const Kokkos::View<DataType, Properties...>& view,
          const ComparatorType& comparator) {...}

template <class ComparatorType, class DataType, class... Properties>
void sort(const Kokkos::View<DataType, Properties...>& view,
          const ComparatorType& comparator) {...}

template <class ExecutionSpace, class ViewType>
std::enable_if_t<Kokkos::is_execution_space<ExecutionSpace>::value> sort(
    const ExecutionSpace& exec, ViewType view, size_t const begin,
    size_t const end) {...}

template <class ViewType>
void sort(ViewType view, size_t const begin, size_t const end) {...}