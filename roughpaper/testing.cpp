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
#include "../libs/sdmgridboxes/movement_in_domain.hpp"
#include "../libs/sdmgridboxes/superdropwithgbxindex.hpp"
#include "../libs/superdrop_solver/sdmmotion.hpp"

namespace dlc = dimless_constants;

void print_nbourmaps(const Maps4GridBoxes &mdlmaps, const double COORD0);
void print_gridboxmaps(const Maps4GridBoxes &mdlmaps, const double COORD0);
void print_superdropcoords(const std::vector<GridBox> &gridboxes,
                           const Maps4GridBoxes &mdlmaps);
int visualise_tsteps();

inline int exchange_or_outstep(const int t_out, const int outstep,
                               const int xchangestep)
/* given current time, t_out, work out which event (exchange or output)
is next to occur and return the time of the sooner event */
{
  const int next_xchange = ((t_out / xchangestep) + 1) * xchangestep; // t of next xchange
  const int next_out = ((t_out / outstep) + 1) * outstep;             // t of next output

  return std::min(next_xchange, next_out);
}

inline int next_step(int t, int step)
{
  return ((t / step) + 1) * step;
}

int main()
{
  const std::string abspath("/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/");
  
  const std::string configfilepath = abspath+"src/config/config.txt";    // path to configuration (.txt file)
  const std::string constantsfilepath = abspath+"src/include/claras_SDconstants.hpp"; // path to constants (.hpp file)
  const Config config(configfilepath, constantsfilepath);

  const std::string grid_filename = abspath+"build/share/dimlessGBxboundaries.dat";    
  const std::string initSDs_filename = abspath+"build/share/dimlessSDsinit.dat";   

  const Maps4GridBoxes mdlmaps(config.SDnspace, grid_filename); 

  const auto solute(std::make_shared<const SoluteProperties>());
  std::vector<SuperdropWithGbxindex>
      SDsInGBxs = create_superdrops_from_initSDsfile(initSDs_filename,
                                              config.nSDsvec,
                                              config.SDnspace, solute);

  /* vector containing all gridboxes that makeup the SDM domain */
  std::vector<GridBox> gridboxes = create_gridboxes(mdlmaps, SDsInGBxs);

  print_gridboxmaps(mdlmaps, dlc::COORD0);
  print_nbourmaps(mdlmaps, dlc::COORD0);
  print_superdropcoords(gridboxes, mdlmaps);

  const double seditstep = 1.0;//timestep2dimlesstime(mdlsteps.sedistep);
  const SdmMotion auto sdmmotion = MoveWithSedimentation(seditstep,
                                                         SimmelTerminalVelocity{});

  move_superdrops_in_domain(mdlmaps, sdmmotion,
                            SDsInGBxs, gridboxes);
  print_superdropcoords(gridboxes, mdlmaps);

  visualise_tsteps();

  return 0;
}
void print_nbourmaps(const Maps4GridBoxes &mdlmaps, const double COORD0)
{
  std::cout << "---- NBOUR MAPS ----\n";
  std::cout << "Z nghbours" << "\n";
  for (const auto idxkey : mdlmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << mdlmaps.get_neighbour_zdown(idxkey) << ", "
              << mdlmaps.get_neighbour_zup(idxkey) << '\n';
  }

  std::cout << "X nghbours" << "\n";
  for (const auto idxkey : mdlmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << mdlmaps.get_neighbour_xbehind(idxkey) << ", "
              << mdlmaps.get_neighbour_xinfront(idxkey) << '\n'; 
  }

  std::cout << "Y nghbours" << "\n";
  for (const auto idxkey : mdlmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << mdlmaps.get_neighbour_yleft(idxkey) << ", "
              << mdlmaps.get_neighbour_yright(idxkey) << '\n';  
  }
  std::cout << "------------------\n";
}

void print_gridboxmaps(const Maps4GridBoxes &mdlmaps, const double COORD0)
{
  std::cout << "---- GBX MAPS ----\n";
  std::cout << "Zmap" << "\n";
  for (const auto idxkey : mdlmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << mdlmaps.get_bounds_z(idxkey).first << ", "
              <<  mdlmaps.get_bounds_z(idxkey).second << '\n'; 
  }

  std::cout << "Xmap" << "\n";
  for (const auto idxkey : mdlmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << mdlmaps.get_bounds_x(idxkey).first << ", "
              <<  mdlmaps.get_bounds_x(idxkey).second << '\n'; 
  }

  std::cout << "Ymap" << "\n";
  for (const auto idxkey : mdlmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << mdlmaps.get_bounds_y(idxkey).first << ", "
              <<  mdlmaps.get_bounds_y(idxkey).second << '\n'; 
  }

  std::cout << "Vol map" << "\n";
  for (const auto idxkey : mdlmaps.gbxidxs)
  {
    std::cout << idxkey << ": "
              << mdlmaps.get_volume(idxkey) 
              << " -> ie. = " << mdlmaps.get_volume(idxkey) * pow(COORD0, 3.0) << "m^3\n";
  }
  std::cout << "----------------\n";
}

void print_superdropcoords(const std::vector<GridBox> &gridboxes,
                           const Maps4GridBoxes &mdlmaps)
{
  std::cout << "\n---- SD Positions -----\n";
  std::cout << " -- in Z direction --\n";
  for (auto gbx: gridboxes)
  {
    std::cout << "GBx " << gbx.gbxindex << " : ("
              << mdlmaps.get_bounds_z(gbx.gbxindex).first << ", " << ", "
              << mdlmaps.get_bounds_z(gbx.gbxindex).second << ")\n";
    
    for (auto SDinGBx : gbx.span4SDsinGBx)
    {
      std::cout << "sdgbx " << SDinGBx.sd_gbxindex
                << " : " << SDinGBx.superdrop.coord3 << "\n";
    }
  }

  std::cout << " -- in X direction --\n";
  for (auto gbx: gridboxes)
  {
    std::cout << "GBx " << gbx.gbxindex << " : ("
              << mdlmaps.get_bounds_x(gbx.gbxindex).first << ", " << ", "
              << mdlmaps.get_bounds_x(gbx.gbxindex).second << ")\n";
    
    for (auto SDinGBx : gbx.span4SDsinGBx)
    {
      std::cout << "sdgbx " << SDinGBx.sd_gbxindex
      << " : " << SDinGBx.superdrop.coord1 << "\n";
    }
  }

  std::cout << " -- Summary --\n";
  for (auto gbx: gridboxes)
  {
    std::cout << "GBx" << gbx.gbxindex << ", ("
              << mdlmaps.get_bounds_x(gbx.gbxindex).first << " " << ", "
              << mdlmaps.get_bounds_x(gbx.gbxindex).second << ") SDs: ";
    
    for (auto SDinGBx : gbx.span4SDsinGBx)
    {
      std::cout << SDinGBx.superdrop.id << ", ";
    }
    std::cout << "\n";
  }

  std::cout << "\n-----------------------\n";
}

int visualise_tsteps()
{
  const int t_end = 10;
  const int couplstep = 3; //outstep
  const int motionstep = 4;
  const int sdmstep = 4;
  int t_out = 0; // time that is incremented by length 'outstep' between coupling communication
  
  while (t_out <= t_end)
  {
    std::cout << "run_sdmstep\n";
    int t = t_out;
    const int nextt = exchange_or_outstep(t, couplstep, motionstep);

    for (int subt = t; subt < nextt; subt=next_step(subt, sdmstep))
    {
      std::cout << "---> process step\n";
    }

    if (nextt % motionstep == 0)
    {
      std::cout << "--> motion step\n";
    }

    std::cout << "run_driverstep\n";
    std::cout << "-> couple step\n";


    t_out += couplstep;
  }


  return 0;
}