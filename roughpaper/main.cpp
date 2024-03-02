/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Monday 29th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 2nd March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * rough paper for checking small things
 */


#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_NestedSort.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "superdrops/collisionkinetics.hpp"

// ---------------------------------------------------------------- //

// struct Superdrop {
//  private:
//   unsigned int sdgbxindex;

//  public:
//   KOKKOS_INLINE_FUNCTION
//   void set_sdgbxindex(const unsigned int ii) { sdgbxindex = ii; }

//   KOKKOS_INLINE_FUNCTION
//   unsigned int get_sdgbxindex() const { return sdgbxindex; }
// };

// using viewd_supers = Kokkos::View<Superdrop *>;
// using viewd_constsupers =
//     Kokkos::View<const Superdrop *>;  // view in device memory of const superdroplets
// using ExecSpace = Kokkos::DefaultExecutionSpace;

// viewd_supers init_supers(const size_t nsupers, const size_t ngbxs) {
//   viewd_supers supers("supers", nsupers);
//   auto h_supers = Kokkos::create_mirror_view(supers);
//   for (size_t kk(0); kk < nsupers; ++kk) {
//     const auto ii = static_cast<unsigned int>(1);
//     h_supers(kk) = Superdrop();
//     h_supers(kk).set_sdgbxindex(ii);
//     std::cout << "ii: " << h_supers(kk).get_sdgbxindex() << "\n";
//   }
//   Kokkos::deep_copy(supers, h_supers);

//   return supers;
// }

// ---------------------------------------------------------------- //

int main(int argc, char *argv[]) {
  const auto r1 = 2.3;
  const auto r2 = 2.3;

  Kokkos::initialize(argc, argv);
  {
    Kokkos::parallel_for(
        "test", 5, KOKKOS_LAMBDA(const unsigned int i) {
          const auto a = coal_surfenergy(r1, r2);
          assert((a == 7.68216e-12) && "a check");
        });
  }
  Kokkos::finalize();
}
