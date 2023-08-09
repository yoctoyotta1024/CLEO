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

void cke_check(const double d1cm,
               const double d2cm,
               const Superdrop &drop1,
               const Superdrop &drop2)
{
  const auto terminalv(SimmelTerminalVelocity{});
  CollisionKinetics ck(terminalv);

  std::cout << "cke check: " << ck.collision_kinetic_energy(drop1, drop2) << "\n";
  
  const double d1(d1cm/100); // [m]
  const double d2(d2cm/100);
  const double vdiff(terminalv(drop1) - terminalv(drop2));
  const double ratio(std::pow(d1, 3.0) * std::pow(d2, 3.0) / (std::pow(d1, 3.0) + std::pow(d2, 3.0)));
  const double ckeconst(dlc::Rho_l * dlc::RHO0 * std::numbers::pi / 12.0);
  const double cke(vdiff * vdiff * ratio * ckeconst);
  std::cout << "breakdown: " << dlc::Rho_l * dlc::RHO0 << ", " << vdiff << ", " << ratio << " -> " << cke << "\n";
}

void surfe_check(const double d1cm,
               const double d2cm,
               const Superdrop &drop1,
               const Superdrop &drop2)
{
  const auto terminalv(SimmelTerminalVelocity{});
  CollisionKinetics ck(terminalv);

  std::cout << "surfe1: " << ck.surfenergy(drop1) << "\n";
  std::cout << "surfe2: " << ck.surfenergy(drop2) << "\n";

  const double d1(d1cm/100); // [m]
  const double d2(d2cm/100);
  const double sigma(7.28e-2);

  std::cout << "1: " << sigma * std::numbers::pi * d1 * d1 << "\n"; 
  std::cout << "2: " << sigma * std::numbers::pi * d2 * d2 << "\n"; 
}

void totsurfe_check(const double d1cm,
               const double d2cm,
               const Superdrop &drop1,
               const Superdrop &drop2)
{
  const auto terminalv(SimmelTerminalVelocity{});
  CollisionKinetics ck(terminalv);

  std::cout << "totsurfe: " << ck.total_surfenergy(drop1, drop2) << "\n";
  
  const double d1(d1cm/100); // [m]
  const double d2(d2cm/100);
  const double sigma(7.28e-2);

  std::cout << "tot: " << sigma * std::numbers::pi * (d1 * d1 + d2 * d2)<< "\n"; 
}

void coalsurfe_check(const double d1cm,
               const double d2cm,
               const Superdrop &drop1,
               const Superdrop &drop2)
{
  const auto terminalv(SimmelTerminalVelocity{});
  CollisionKinetics ck(terminalv);

  std::cout << "totsurfe: " << ck.coal_surfenergy(drop1, drop2) << "\n";
  
  const double d1(d1cm/100); // [m]
  const double d2(d2cm/100);
  const double sigma(7.28e-2);

  const double d3sum = std::pow((d1 * d1 * d1 + d2 * d2 * d2), 2.0/3.0);
  std::cout << "tot: " << sigma * std::numbers::pi * d3sum << "\n"; 
}

int main()
{
  const double d1cm(0.18); // [cm]
  const double d2cm(0.1); // [cm]

  auto sdIdGen = Superdrop::IDType::Gen{};
  
  const auto isolute(std::make_shared<const SoluteProperties>());
  const double im_sol = 3.0;

  const unsigned long long ieps1 = 3;
  const double iradius1 = d1cm / 200.0 / dlc::R0;
  Superdrop drop1(isolute, ieps1, iradius1, im_sol, 0.0, 0.0, 0.0, sdIdGen.next());

  const unsigned long long ieps2 = 6;
  const double iradius2 = d2cm / 200.0 / dlc::R0;;
  Superdrop drop2(isolute, ieps2, iradius2, im_sol, 0.0, 0.0, 0.0, sdIdGen.next());

  cke_check(d1cm, d2cm, drop1, drop2);
  surfe_check(d1cm, d2cm, drop1, drop2);
  totsurfe_check(d1cm, d2cm, drop1, drop2);
  coalsurfe_check(d1cm, d2cm, drop1, drop2);

  auto ll = LowListCollCoalEff(SimmelTerminalVelocity{});
  std::cout <<"\ncoaleff: " << ll.coaleff(drop1, drop2) << "\n";

  return 0;
}