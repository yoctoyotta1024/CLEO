// Author: Clara Bayley
// File: thermodynamicsfromfile.cpp
/* implementation of things specifically
to run uncoupled SDM where thermodynamics
are read from file */

#include "thermodynamicsfromfile.hpp"

std::vector<double>
thermodynamicvar_from_binary(std::string_view filename)
{
  /* open file and read in the metatdata
  for all the variables in that file */
  std::ifstream file(open_binary(filename));
  std::vector<VarMetadata> meta(metadata_from_binary(file));
  
  /* read in the data for the 1st variable in the file */
  std::vector<double>
      thermovar(vector_from_binary<double>(file, meta.at(0)));

  return thermovar;
}

std::array<size_t, 3> kijfromindex(const std::array<size_t, 3> &ndims,
                                   const size_t index)
/* return (k,i,j) indicies from idx for a flattened 3D array
with ndims [nz, nx, ny]. kij is useful for then getting
position in of a variable in a flattened array defined on
the faces of the same grid. E.g for the w velocity defined 
on z faces of the grid which therefore has dims [nz+1, nx, ny] */
{
  const size_t j = index / (ndims[0] * ndims[1]);
  const size_t k = index % ndims[0];
  const size_t i = index / ndims[0] - ndims[1] * j;

  return std::array<size_t, 3>{k, i, j};
}

ThermodynamicsFromFile::
    ThermodynamicsFromFile(const Config &config,
                           const std::array<size_t, 3> &ndims,
                           const size_t nsteps)
    : ndims(ndims),
      atpos(0),
      atpos_zface(0),
      atpos_xface(0),
      atpos_yface(0),
      press(thermodynamicvar_from_binary(config.press_filename)),
      temp(thermodynamicvar_from_binary(config.temp_filename)),
      qvap(thermodynamicvar_from_binary(config.qvap_filename)),
      qcond(thermodynamicvar_from_binary(config.qcond_filename)),
      wvel(), uvel(), vvel(),
      get_wvelzfaces([](const unsigned int ii)
                     { return std::pair<double, double>{0.0, 0.0}; }),
      get_uvelxfaces([](const unsigned int ii)
                     { return std::pair<double, double>{0.0, 0.0}; }),
      get_vvelyfaces([](const unsigned int ii)
                     { return std::pair<double, double>{0.0, 0.0}; })
{
  std::string windstr = set_windvelocities(config);
  std::cout << "\nFinished reading thermodynamics from binaries for:\n"
               "  pressure,\n  temperature,\n"
               "  water vapour mass mixing ratio,\n"
               "  liquid water mass mixing ratio,\n  "
            << windstr << '\n';

  check_thermodyanmics_vectorsizes(config.SDnspace, ndims, nsteps);
}

std::string ThermodynamicsFromFile::
    set_windvelocities(const Config &config)
/* depending on SDnspace, read in data
for wind velocity components (or not)*/
{
  const int SDnspace(config.SDnspace);

  if (SDnspace == 0)
  {
    return "0-D model has no wind data";
  }
  else if (SDnspace <= 3) // means 1 <= SDnspace < 4
  {
    return set_windvelocities_frombinaries(config);
  }
  else // means SDnspace > 3
  {
    const std::string errmsg("SDnspace > 3 is invalid");
    throw std::invalid_argument(errmsg);
  }
}

std::string ThermodynamicsFromFile::
    set_windvelocities_frombinaries(const Config &config)
/* Read in data from binary files for wind
velocity components in 1D, 2D or 3D model
and check they have correct size */
{
  const int SDnspace(config.SDnspace);
 
  std::string info(std::to_string(SDnspace) + "-D model ");
  if (SDnspace >= 1)
  {
    wvel = wvel_from_binary(config.wvel_filename);
    if (SDnspace >= 2)
    {
      uvel = uvel_from_binary(config.uvel_filename);
      if (SDnspace == 3)
      {
        vvel = vvel_from_binary(config.vvel_filename);
        return info + "[w, u, v] wind velocity";
      }
      return info+"[w, u] wind velocity";
    }
    return info+"vertical w wind velocity";
  }
  return info;
}

std::vector<double> ThermodynamicsFromFile::
    wvel_from_binary(std::string_view filename)
/* set function for retrieving wvel defined at zfaces of
a gridbox with index 'gbxindex' and return vector
containting wvel data from binary file */
{
  get_wvelzfaces = [&](const unsigned int gbxindex)
  {
    const auto kij = kijfromindex(ndims, (size_t)gbxindex); // [k,i,j] of gridbox centre on 3D grid
    const size_t nzfaces(ndims[0] + 1);                     // no. z faces to same 3D grid

    size_t lpos(ndims[1] * nzfaces * kij[2] + nzfaces * kij[1] + kij[0]); // position of z lower face in 1D wvel vector
    lpos += atpos_zface;
    const size_t uppos(lpos + 1); // position of z upper face

    return std::pair<double, double>{wvel.at(lpos), wvel.at(uppos)};
  };

  return thermodynamicvar_from_binary(filename);
}

std::vector<double> ThermodynamicsFromFile::
    uvel_from_binary(std::string_view filename)
/* set function for retrieving uvel defined at xfaces of
a gridbox with index 'gbxindex' and return vector
containting uvel data from binary file */
{
  get_uvelxfaces = [&](const unsigned int gbxindex)
  {
    const auto kij = kijfromindex(ndims, (size_t)gbxindex); // [k,i,j] of gridbox centre on 3D grid
    const size_t nxfaces(ndims[1] + 1);                     // no. x faces to same 3D grid

    size_t lpos(nxfaces * ndims[0] * kij[2] + ndims[0] * kij[1] + kij[0]); // position of x lower face in 1D uvel vector
    lpos += atpos_xface;
    const size_t uppos(lpos + ndims[0]); // position of x upper face

    return std::pair<double, double>{uvel.at(lpos), uvel.at(uppos)};
  };

  return thermodynamicvar_from_binary(filename);
}

std::vector<double> ThermodynamicsFromFile::
    vvel_from_binary(std::string_view filename)
    /* set function for retrieving vvel defined at yfaces of 
    a gridbox with index 'gbxindex' and return vector 
    containting vvel data from binary file */
{
  get_vvelyfaces = [&](const unsigned int gbxindex)
  {
    const auto kij = kijfromindex(ndims, (size_t)gbxindex); // [k,i,j] of gridbox centre on 3D grid
    
    size_t lpos(ndims[1] * ndims[0] * kij[2] + ndims[0] * kij[1] + kij[0]); // position of y lower face in 1D vvel vector
    lpos += atpos_yface;
    const size_t uppos(lpos + ndims[1] * ndims[0]); // position of x upper face

    return std::pair<double, double>{vvel.at(lpos), vvel.at(uppos)};
  };

  return thermodynamicvar_from_binary(filename);
}

void ThermodynamicsFromFile::
    check_thermodyanmics_vectorsizes(const int SDnspace,
                                     const std::array<size_t, 3> &ndims,
                                     const size_t nsteps) const
{
  auto is_size = [](const std::vector<double> &vel, const size_t sz)
  { 
    const size_t velsize(vel.size());
    if( velsize != sz )
    {
      throw std::invalid_argument(std::to_string(velsize)+ " vector is "
                                "not consistent with correct size "+
                                std::to_string(sz));
    }
  };

  const size_t ngridboxes(ndims[0]*ndims[1]*ndims[2]);
  is_size(press, nsteps*ngridboxes);
  check_vectorsizes({press.size(), temp.size(),
                      qvap.size(), qcond.size()});

  if (SDnspace >= 1)
  {
    const size_t wsz = nsteps*(ndims[0]+1)*ndims[1]*ndims[2];
    is_size(wvel, wsz);
    if (SDnspace >= 2)
    {
      const size_t usz = nsteps*ndims[0]*(ndims[1]+1)*ndims[2];
      is_size(uvel, usz);
      if (SDnspace >= 3)
      {
        const size_t vsz = nsteps*ndims[0]*ndims[1]*(ndims[2]+1);
        is_size(vvel, vsz);
      }
    } 
  }
}