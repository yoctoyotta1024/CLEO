/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Monday 29th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 6th March 2024
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

#include "src/collisionkinetics.hpp"

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
