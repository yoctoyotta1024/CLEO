// Author: Clara Bayley
// File: superdropwithgbxindex.cpp
/* functionality involved in handling vector of
SuperdropWithGbxindex instances (see superdrop.hpp
for definition of this struct).
Note: some functions are declared here to avoid
being visible externally */

#include "superdropwithgbxindex.hpp"

/* -- function called by create_superdrops_from_initSDsfile -- */
Kokkos::vector<SuperdropWithGbxindex>
create_superdropwithgbxindexes(const int nSDsvec, const int SDnspace,
                               const InitSDsData &initSDs,
                               const std::shared_ptr<const SoluteProperties>
                                   solute);

std::vector<double> initSDcoords(const int SDnspace,
                                 const InitSDsData &initSDs,
                                 const int i);
/* ------------------------------------------------------ */

Kokkos::vector<SuperdropWithGbxindex>
create_superdrops_from_initSDsfile(std::string_view initSDs_filename,
                            const int nSDsvec,
                            const int SDnspace,
                            const std::shared_ptr<const SoluteProperties> solute)
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

  Kokkos::vector<SuperdropWithGbxindex>
      SDsInGBxs = create_superdropwithgbxindexes(nSDsvec, SDnspace,
                                                 initSDs, solute);

  /* 3. Initialise gridbox index associated with each superdroplets */
  std::cout << "Now sorting superdroplets based on the"
               " index of their associated gridboxes\n";

  sort_superdrops_via_gridboxindex(SDsInGBxs);

  return SDsInGBxs;
}

Kokkos::vector<SuperdropWithGbxindex>
create_superdropwithgbxindexes(const int nSDsvec, const int SDnspace,
                               const InitSDsData &initSDs,
                               const std::shared_ptr<const SoluteProperties>
                                   solute)
{
  Kokkos::vector<SuperdropWithGbxindex> SDsInGBxs;
  auto sdIdGen = Superdrop::IDType::Gen{};

  for (int i = 0; i < nSDsvec; ++i)
  {
    const unsigned int sd_gbxindex = initSDs.sd_gbxindex.at(i);
    const auto sd_identity = sdIdGen.next();
    const unsigned long long eps = (unsigned long long)(initSDs.eps_init.at(i) + 0.5);
    const double radius = initSDs.radius_init.at(i);
    const double m_sol = initSDs.m_sol_init.at(i);
    const std::vector<double> zxycoords = initSDcoords(SDnspace, initSDs, i);

    const SuperdropWithGbxindex SDinGBx(sd_gbxindex,
                                       Superdrop(solute, eps, radius, m_sol,
                                                 zxycoords.at(0), zxycoords.at(1),
                                                 zxycoords.at(2), sd_identity));
    
    // print_SDinGBx(SDinGBx);
    SDsInGBxs.push_back(SDinGBx);
  }

  if (SDsInGBxs.size() < initSDs.sd_gbxindex.size())
  {
    const std::string err("Fewer superdroplets were created than were"
                          " read from initialisation file into initSDs"
                          " ie. " +
                          std::to_string(SDsInGBxs.size()) + " < " +
                          std::to_string(initSDs.sd_gbxindex.size()));
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