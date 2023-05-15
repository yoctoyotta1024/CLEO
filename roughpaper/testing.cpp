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
#include <vector>
#include <cmath>
#include <iostream>
#include <string_view>

#include <Kokkos_core.hpp>

#include "../libs/claras_SDconstants.hpp"
#include "../libs/initialisation/config.hpp"
#include "../libs/initialisation/readbinary.hpp"
#include "../libs/initialisation/read_gbxboundaries.hpp"
#include "../libs/initialisation/read_initsuperdrops.hpp"
#include "../libs/sdmgridboxes/maps4gridboxes.hpp"
#include "../libs/sdmgridboxes/gridbox.hpp"
#include "../libs/sdmgridboxes/movesuperdropsindomain.hpp"
#include "../libs/sdmgridboxes/superdropwithgbxindex.hpp"
#include "../libs/sdmgridboxes/runsdmstep.hpp"
#include "../libs/sdmgridboxes/sdmotion.hpp"

namespace dlc = dimless_constants;

void print_nbourmaps(const Maps4GridBoxes &gbxmaps, const double COORD0);
void print_gridboxmaps(const Maps4GridBoxes &gbxmaps, const double COORD0);
void print_superdropcoords(const std::vector<GridBox> &gridboxes,
                           const Maps4GridBoxes &gbxmaps);

int main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  Kokkos::DefaultExecutionSpace{}.print_configuration(std::cout);

  const long n = 101;

  std::cout << "Number of even integers from 0 to " << n-1 << "\n";

  Kokkos::Timer timer;
  timer.reset();

  // Compute the number of even integers from 0 to n-1, in parallel.
  long count = 0;
  Kokkos::parallel_reduce(
      n, KOKKOS_LAMBDA(const long i, long& lcount) { lcount += (i % 2) == 0; },
      count);

  double count_time = timer.seconds();
  std::cout << "  Parallel: "<<count << "    " << count_time << "\n";

  timer.reset();

  // Compare to a sequential loop.
  long seq_count = 0;
  for (long i = 0; i < n; ++i) {
    seq_count += (i % 2) == 0;
  }

  count_time = timer.seconds();
  std::cout << "  Sequential: "<<seq_count << "    " << count_time << "\n";

  Kokkos::finalize();

  return (count == seq_count) ? 0 : -1;
}


void print_nbourmaps(const Maps4GridBoxes &gbxmaps, const double COORD0)
{
  std::cout << "---- NBOUR MAPS ----\n";
  std::cout << "Z nghbours" << "\n";
  for (const auto idxkey : gbxmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << gbxmaps.get_neighbour_zdown(idxkey) << ", "
              << gbxmaps.get_neighbour_zup(idxkey) << '\n';
  }

  std::cout << "X nghbours" << "\n";
  for (const auto idxkey : gbxmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << gbxmaps.get_neighbour_xbehind(idxkey) << ", "
              << gbxmaps.get_neighbour_xinfront(idxkey) << '\n'; 
  }

  std::cout << "Y nghbours" << "\n";
  for (const auto idxkey : gbxmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << gbxmaps.get_neighbour_yleft(idxkey) << ", "
              << gbxmaps.get_neighbour_yright(idxkey) << '\n';  
  }
  std::cout << "------------------\n";
}

void print_gridboxmaps(const Maps4GridBoxes &gbxmaps, const double COORD0)
{
  std::cout << "---- GBX MAPS ----\n";
  std::cout << "Zmap" << "\n";
  for (const auto idxkey : gbxmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << gbxmaps.get_bounds_z(idxkey).first << ", "
              <<  gbxmaps.get_bounds_z(idxkey).second << '\n'; 
  }

  std::cout << "Xmap" << "\n";
  for (const auto idxkey : gbxmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << gbxmaps.get_bounds_x(idxkey).first << ", "
              <<  gbxmaps.get_bounds_x(idxkey).second << '\n'; 
  }

  std::cout << "Ymap" << "\n";
  for (const auto idxkey : gbxmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << gbxmaps.get_bounds_y(idxkey).first << ", "
              <<  gbxmaps.get_bounds_y(idxkey).second << '\n'; 
  }

  std::cout << "Vol map" << "\n";
  for (const auto idxkey : gbxmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << gbxmaps.get_volume(idxkey) 
              << " -> ie. = " << gbxmaps.get_volume(idxkey) * pow(COORD0, 3.0) << "m^3\n";
  }
  std::cout << "----------------\n";
}

void print_superdropcoords(const std::vector<GridBox> &gridboxes,
                           const Maps4GridBoxes &gbxmaps)
{
  std::cout << "\n---- SD Positions -----\n";
  std::cout << " -- in Z direction --\n";
  for (auto gbx: gridboxes)
  {
    std::cout << "GBx " << gbx.gbxindex << " : ("
              << gbxmaps.get_bounds_z(gbx.gbxindex).first << ", " << ", "
              << gbxmaps.get_bounds_z(gbx.gbxindex).second << ")\n";
    
    for (auto SDinGBx : gbx.span4SDsinGBx)
    {
      std::cout << "sdgbx " << SDinGBx.sd_gbxindex
                << " : " << SDinGBx.superdrop.coord3 << "\n";
    }
  }

  std::cout << " -- Summary --\n";
  for (auto gbx: gridboxes)
  {
    std::cout << "GBx" << gbx.gbxindex << ", ("
              << gbxmaps.get_bounds_x(gbx.gbxindex).first << " " << ", "
              << gbxmaps.get_bounds_x(gbx.gbxindex).second << ") SDs: ";
    
    for (auto SDinGBx : gbx.span4SDsinGBx)
    {
      std::cout << SDinGBx.superdrop.id << ", ";
    }
    std::cout << "\n";
  }

  std::cout << "\n-----------------------\n";
}