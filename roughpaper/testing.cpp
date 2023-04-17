// Author: Clara Bayley
// File: testing.cpp
/* This file runs test snippets of c++.
could compile with e.g.
/opt/homebrew/bin/g++-12 testing.cpp ../src/include/superdrop_solver/superdrop.cpp -I ../src/include/superdrop_solver --std=c++20
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

#include "../libs/claras_SDconstants.hpp"
#include "../libs/initialisation/config.hpp"
#include "../libs/initialisation/readbinary.hpp"
#include "../libs/initialisation/read_gbxboundaries.hpp"
#include "../libs/initialisation/read_initsuperdrops.hpp"
#include "../libs/sdmgridboxes/maps4gridboxes.hpp"
#include "../libs/sdmgridboxes/gridbox.hpp"
#include "../libs/sdmgridboxes/movesuperdropsindomain.hpp"
#include "../libs/sdmgridboxes/superdropwithgbxindex.hpp"
#include "../libs/sdmgridboxes/run_sdmstep.hpp"
#include "../libs/superdrop_solver/sdmotion.hpp"

namespace dlc = dimless_constants;

void print_nbourmaps(const Maps4GridBoxes &gbxmaps, const double COORD0);
void print_gridboxmaps(const Maps4GridBoxes &gbxmaps, const double COORD0);
void print_superdropcoords(const std::vector<GridBox> &gridboxes,
                           const Maps4GridBoxes &gbxmaps);

int main()
{
  const std::string abspath("/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/");
  
  const std::string configfilepath = abspath+"src/config/config.txt";    // path to configuration (.txt file)
  const std::string constantsfilepath = abspath+"src/include/claras_SDconstants.hpp"; // path to constants (.hpp file)
  const Config config(configfilepath, constantsfilepath);

  const std::string grid_filename = abspath+"build/share/dimlessGBxboundaries.dat";    
  const std::string initSDs_filename = abspath+"build/share/dimlessSDsinit.dat";   

  const Maps4GridBoxes gbxmaps(config.SDnspace, grid_filename); 

  const auto solute(std::make_shared<const SoluteProperties>());
  std::vector<SuperdropWithGbxindex>
      SDsInGBxs = create_superdrops_from_initSDsfile(initSDs_filename,
                                              config.nSDsvec,
                                              config.SDnspace, solute);

  /* vector containing all gridboxes that makeup the SDM domain */
  std::vector<GridBox> gridboxes = create_gridboxes(gbxmaps, SDsInGBxs);

  print_gridboxmaps(gbxmaps, dlc::COORD0);
  print_nbourmaps(gbxmaps, dlc::COORD0);
  print_superdropcoords(gridboxes, gbxmaps);

  const SdMotion auto movesd = NullMotion();
  const MoveSuperdropsInDomain move(movesd);

  move.move_superdrops_in_domain(gbxmaps, SDsInGBxs, gridboxes);
  print_superdropcoords(gridboxes, gbxmaps);

  return 0;
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