// Author: Clara Bayley
// File: testing.cpp
/* This file runs test snippets of c++.
could compile with e.g.
/opt/homebrew/bin/g++-13 testing.cpp ../src/include/superdrop_solver/superdrop.cpp -I ../src/include/superdrop_solver --std=c++20
*/

#include <map>
#include <limits>
#include <ranges>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string_view>
#include <span>
#include <utility> 

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>

#include "../libs/claras_SDconstants.hpp"
// #include "../libs/initialisation/config.hpp"
// #include "../libs/initialisation/readbinary.hpp"
// #include "../libs/initialisation/read_gbxboundaries.hpp"
// #include "../libs/initialisation/read_initsuperdrops.hpp"
// #include "../libs/sdmgridboxes/maps4gridboxes.hpp"
// #include "../libs/sdmgridboxes/gridbox.hpp"
// #include "../libs/sdmgridboxes/movesuperdropsindomain.hpp"
// #include "../libs/sdmgridboxes/superdropwithgbxindex.hpp"
// #include "../libs/sdmgridboxes/runsdmstep.hpp"
// #include "../libs/sdmgridboxes/sdmotion.hpp"
#include "../libs/superdrop_solver/collisionxkernels.hpp"

namespace dlc = dimless_constants;

int main()
{
  auto sdIdGen = Superdrop::IDType::Gen{};
  
  const auto isolute(std::make_shared<const SoluteProperties>());
  const double im_sol = 3.0;
  
  const unsigned long long ieps1 = 3;
  const double iradius1 = 1.0;
  Superdrop drop1(isolute, ieps1, iradius1, im_sol, 0.0, 0.0, 0.0, sdIdGen.next());

  const unsigned long long ieps2 = 6;
  const double iradius2 = 0.25;
  Superdrop drop2(isolute, ieps2, iradius2, im_sol, 0.0, 0.0, 0.0, sdIdGen.next());

  std::cout <<"drop1: " << drop1.radius <<"\ndrop2: " << drop2.radius << "\n";

  const auto llke(LowListKernelEfficiency{SimmelTerminalVelocity{}});

  std::cout << "llke: " << llke(drop1, drop2) <<"\n";
  return 0;
}