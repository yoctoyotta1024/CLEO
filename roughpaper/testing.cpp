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

void print_gridboxmaps(const Maps4GridBoxes &mdlmaps, const double COORD0);
void print_initSDs(const InitSDsData &initSDs);

int main()
{
  const std::string abspath("/Users/yoctoyotta1024/Documents/autumnwinter2022_23"
                            "/clara-bayley-superdroplets/superdroplet_model/");
  
  const std::string configfilepath = abspath+"src/config/config.txt";    // path to configuration (.txt file)
  const std::string constantsfilepath = abspath+"src/include/claras_SDconstants.hpp"; // path to constants (.hpp file)
  const Config config(configfilepath, constantsfilepath);

  const std::string grid_filename = abspath+"build/"+config.grid_filename;
  const Maps4GridBoxes mdlmaps(config.SDnspace, grid_filename); 

  print_gridboxmaps(mdlmaps, dlc::COORD0);

  const std::string initSDs_filename = abspath+"build/"+config.initSDs_filename;
  const InitSDsData initSDs = get_initsuperdropsdata(initSDs_filename);

  print_initSDs(initSDs);

  const auto solute(std::make_shared<const SoluteProperties>());
  std::vector<SuperdropWithGridbox>
      SDsInGBxs = superdrops_from_initSDsfile(initSDs_filename,
                                              config.nSDsvec,
                                              config.SDnspace, solute,
                                              mdlmaps);

  for (auto a : SDsInGBxs)
  {
    std::cout << "---\nSD " << a.superdrop.id.value
              << ": " << a.sd_gbxindex << ", " << a.superdrop.eps
              << ", " << a.superdrop.radius << ", " << a.superdrop.m_sol
              << ", " << a.superdrop.coord3 << ", " << a.superdrop.coord1
              << ", " << a.superdrop.coord2 << "\n ---- ";
  }

  return 0;
}

void print_gridboxmaps(const Maps4GridBoxes &mdlmaps, const double COORD0)
{
  std::cout << "---- RESULTS ----\n";
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

  std::cout << "m_sol_init" << "\n";
  for (auto x : initSDs.m_sol_init)
  {
    std::cout << x << ", ";
  }
  std::cout << "\n----------------\n";

  std::cout << "coord3_init" << "\n";
  for (auto x : initSDs.coord3_init)
  {
    std::cout << x << ", ";
  }
  std::cout << "\n----------------\n";

  std::cout << "coord1_init" << "\n";
  for (auto x : initSDs.coord1_init)
  {
    std::cout << x << ", ";
  }
  std::cout << "\n----------------\n";

  std::cout << "coord2_init" << "\n";
  for (auto x : initSDs.coord2_init)
  {
    std::cout << x << ", ";
  }
  std::cout << "\n----------------\n";
}