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
#include "../src/runmodel/maps4gridboxes.hpp"
#include "../src/runmodel/superdrops_with_gridboxes.hpp"

namespace dlc = dimless_constants;

void print_nbourmaps(const Maps4GridBoxes &mdlmaps, const double COORD0);
void print_gridboxmaps(const Maps4GridBoxes &mdlmaps, const double COORD0);
void print_initSDs(const InitSDsData &initSDs);

int main()
{
  const std::string abspath("/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/");
  
  const std::string configfilepath = abspath+"src/config/config.txt";    // path to configuration (.txt file)
  const std::string constantsfilepath = abspath+"src/include/claras_SDconstants.hpp"; // path to constants (.hpp file)
  const Config config(configfilepath, constantsfilepath);

  const std::string grid_filename = abspath+"build/"+config.grid_filename;
  const Maps4GridBoxes mdlmaps(config.SDnspace, grid_filename); 

  print_gridboxmaps(mdlmaps, dlc::COORD0);

  const std::string initSDs_filename = abspath+"build/"+config.initSDs_filename;
  const InitSDsData initSDs = get_initsuperdropsdata(initSDs_filename);

  // print_initSDs(initSDs);
  print_nbourmaps(mdlmaps, dlc::COORD0);

  return 0;
}
void print_nbourmaps(const Maps4GridBoxes &mdlmaps, const double COORD0)
{
  std::cout << "---- NBOUR MAPS ----\n";
  std::cout << "Z nghbours" << "\n";
  for (const auto & [key, value] : mdlmaps.idx2nghbour_z)
  {
    std::cout << key << ": " << value.first  <<", " << value.second << '\n';
  }

  std::cout << "X nghbours" << "\n";
  for (const auto & [key, value] : mdlmaps.idx2nghbour_x)
  {
    std::cout << key << ": " << value.first  <<", " << value.second << '\n';
  }

  std::cout << "Y nghbours" << "\n";
  for (const auto & [key, value] : mdlmaps.idx2nghbour_y)
  {
    std::cout << key << ": " << value.first  <<", " << value.second << '\n';
  }
}

void print_gridboxmaps(const Maps4GridBoxes &mdlmaps, const double COORD0)
{
  std::cout << "---- GBX MAPS ----\n";
  std::cout << "Zmap" << "\n";
  for (const auto & [key, value] : mdlmaps.idx2bounds_z)
  {
    std::cout << key << ": " << value.first  <<", " << value.second << '\n';
  }

  std::cout << "Xmap" << "\n";
  for (const auto & [key, value] : mdlmaps.idx2bounds_x)
  {
    std::cout << key << ": " << value.first  <<", " << value.second << '\n';
  }

  std::cout << "Ymap" << "\n";
  for (const auto & [key, value] : mdlmaps.idx2bounds_y)
  {
    std::cout << key << ": " << value.first  <<", " << value.second << '\n';
  }

  std::cout << "Vol map" << "\n";
  for (const auto & [key, value] : mdlmaps.idx2vol)
  {
    std::cout << key << ": " << value << " -> ie. = "<< value * pow(COORD0, 3.0) << "m^3\n";
  }
  std::cout << "----------------\n";
}

void print_initSDs(const InitSDsData &initSDs)
{
  std::cout << "eps_init" << "\n";
  for (auto x : initSDs.eps_init)
  {
    std::cout << x << ", ";
  }
  std::cout << "\n----------------\n";

  std::cout << "radius_init" << "\n";
  for (auto x : initSDs.radius_init)
  {
    std::cout << x << ", ";
  }
  std::cout << "\n----------------\n";


}