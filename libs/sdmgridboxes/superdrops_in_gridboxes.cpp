// Author: Clara Bayley
// File: superdrops_in_gridboxes.cpp
/* functionality involved in handling
the SuperdropWithGbxindex instances (see superdrop.hpp
for definition of this struct). Some functions
declared here to avoid being visible externally */

#include "superdrops_in_gridboxes.hpp"

/* ------------------------------------------------------ */
/* -- function called create_superdrops_in_gridboxes -- */
std::vector<SuperdropWithGbxindex>
create_superdropsingridboxes(const int nSDsvec, const int SDnspace,
                               const InitSDsData &initSDs,
                               const std::shared_ptr<const SoluteProperties> solute,
                               const Maps4GridBoxes &mdlmaps);

std::vector<double> initSDcoords(const int SDnspace,
                                 const InitSDsData &initSDs,
                                 const int i);
/* ------------------------------------------------------ */
/* ----- function called by sdgbxindex_to_neighbour ----- */
int flag_tochange_sdgbxindex(const SuperdropWithGbxindex &SDinGBx,
                             const std::map<unsigned int,
                                            std::pair<double, double>> &idx2bounds_z);
/* ------------------------------------------------------ */

std::vector<SuperdropWithGbxindex>
superdrops_from_initSDsfile(std::string_view initSDs_filename,
                            const int nSDsvec,
                            const int SDnspace,
                            const std::shared_ptr<const SoluteProperties> solute,
                            const Maps4GridBoxes &mdlmaps)
/* reads initsuperdrop file for superdroplets' initial properties. Uses this data
to create 'nSDsvec' no. of SuperdropletWithGridbox instances in a vector
where all the superdroplets have the same solute properties, "solute".
Uses the coordinates of each superdroplet to set the value of the sd_gbxindex
associated with each superdroplet in the SuperdropletWithGridbox struct */
{
  /* 1. Read initial superdroplets' data from 'initsuperdrop' file */
  const InitSDsData initSDs = get_initsuperdropsdata(initSDs_filename);

  /* 2. Create vector of 'nSDsvec' elements. Each element is
  superdroplet with the index of its associated gridbox */
  std::cout << "Initialisation data for superdrops' read from "
            << initSDs_filename << ". "
            << "\nNow creating superdrops with gridboxes\n";

  std::vector<SuperdropWithGbxindex>
      SDsInGBxs = create_superdropsingridboxes(nSDsvec, SDnspace, initSDs,
                                                 solute, mdlmaps);

  /* 3. Initialise gridbox index associated with each superdroplets */
  std::cout << "Now sorting superdroplets based on the"
               " index of their associated gridboxes\n";

  sort_superdrops_via_gridboxindex(SDsInGBxs);

  return SDsInGBxs;
}

std::vector<SuperdropWithGbxindex>
create_superdropsingridboxes(const int nSDsvec, const int SDnspace,
                               const InitSDsData &initSDs,
                               const std::shared_ptr<const SoluteProperties> solute,
                               const Maps4GridBoxes &mdlmaps)
{
  std::vector<SuperdropWithGbxindex> SDsInGBxs;
  auto sdIdGen = Superdrop::IDType::Gen{};

  for (int i = 0; i < nSDsvec; ++i)
  {
    const unsigned int sd_gbxindex = initSDs.sd_gbxindex.at(i);
    const auto sd_identity = sdIdGen.next();
    const size_t eps = (size_t)(initSDs.eps_init.at(i) + 0.5);
    const double radius = initSDs.radius_init.at(i);
    const double m_sol = initSDs.m_sol_init.at(i);
    const std::vector<double> zxycoords = initSDcoords(SDnspace, initSDs, i);

    const SuperdropWithGbxindex SDinGBx(sd_gbxindex,
                                       Superdrop(solute, eps, radius, m_sol,
                                                 zxycoords.at(0), zxycoords.at(1),
                                                 zxycoords.at(2), sd_identity));
    
    print_SDinGBx(SDinGBx);
    SDsInGBxs.push_back(SDinGBx);
  }

  if (SDsInGBxs.size() < initSDs.sd_gbxindex.size())
  {
    const std::string err = "Fewer superdroplets were created than were"
                            " read from initialisation file into initSDs"
                            " ie. "+std::to_string(SDsInGBxs.size())+
                            " < "+std::to_string(initSDs.sd_gbxindex.size());
    throw std::invalid_argument(err);
  }

  return SDsInGBxs;
}

std::vector<double> initSDcoords(const int SDnspace,
                                 const InitSDsData &initSDs,
                                 const int i)
{
  std::vector<double> zxycoords = {0.0, 0.0, 0.0};

  if (SDnspace >= 1)
  {
    const double coord3 = initSDs.coord3_init.at(i);
    zxycoords.at(0) = coord3;
  }
  if (SDnspace >= 2)
  {
    const double coord1 = initSDs.coord1_init.at(i);
    zxycoords.at(1) = coord1;
  }
  if (SDnspace == 3)
  {
    const double coord2 = initSDs.coord2_init.at(i);
    zxycoords.at(2) = coord2;
  }

  return zxycoords;
}

int flag_tochange_sdgbxindex(const SuperdropWithGbxindex &SDinGBx,
                             const std::map<unsigned int, std::pair<double, double>> &idx2bounds_z)
/* Determines value of the is_change flag used to signal if
the gridboxindex associated with a superdrop needs to change
and if so, in which direction the superdroplet needs to move */
{
  const double coord = SDinGBx.superdrop.coord3;
  const double llim = (*idx2bounds_z.find(SDinGBx.sd_gbxindex)).second.first;
  const double ulim = (*idx2bounds_z.find(SDinGBx.sd_gbxindex)).second.second;

  if (coord < llim)
  {
    return 1; // enum for moving SD index down a gridbox
  }
  else if (coord >= ulim)
  {
    return 2; // enum for moving SD index up a gridbox
  }
  else
  {
    return 0; // enum for not changing SD index
  }
}

void sdgbxindex_to_neighbour(const Maps4GridBoxes &mdlmaps,
                             SuperdropWithGbxindex &SDinGBx)
/* first check if gridbox index associated with the superdrop
in SDinGBx needs to change. If it does, implement change by
calling correct function for changing the sd_gbxindex to a
neighbouring gridbox's index in a particular direction.
The direction is given by the value of the is_change flag */
{
  enum MoveSuperdrop
  {
    No,
    Down,
    Up,
    Left,
    Right,
    Forwards,
    Backwards
  };

  int is_change = flag_tochange_sdgbxindex(SDinGBx, mdlmaps.idx2bounds_z);
  
  if (is_change != No)
  {
    if (is_change == Down)
    {
      SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_zdown(SDinGBx.sd_gbxindex);
    }
    else if (is_change == Up)
    {
      SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_zup(SDinGBx.sd_gbxindex);
    }
    else
    {
      throw std::invalid_argument("method to change SD index for"
                                  "given is_change flag is not defined");
    }
  }
}