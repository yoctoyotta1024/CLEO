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

ThermodynamicsFromFile::
    ThermodynamicsFromFile(const Config &config,
                           const std::array<size_t, 3> &ndims,
                           const size_t nsteps)
    : ndims(ndims),
      atpos(0),
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
                     { return std::pair<double, double>{0.0, 0.0}; }),
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
/* depending on SDnspace, read in data for wind velocity
components from binary files into appropriate vectors
and check they have correct size */
{
  const int SDnspace(config.SDnspace);

  if (SDnspace == 0)
  {
    return "0-D model has no wind data";
  }

  else if (SDnspace <= 3) // means 1 <= SDnspace < 4
  {
    std::string info(std::to_string(SDnspace) + "-D model ");
    wvel = thermodynamicvar_from_binary(config.wvel_filename);
    get_wvelzfaces = [&](const unsigned int gbxindex)
    { return wvel.at(atpos + (size_t)gbxindex); };

    if (SDnspace >= 2)
    {
      uvel = thermodynamicvar_from_binary(config.uvel_filename);
      get_uvelxfaces = [&](const unsigned int gbxindex)
      { return uvel.at(atpos + (size_t)gbxindex); };

      if (SDnspace == 3)
      {
        vvel = thermodynamicvar_from_binary(config.vvel_filename);
        get_vvelyfaces = [&](const unsigned int gbxindex)
        { return vvel.at(atpos + (size_t)gbxindex); };
        
        return info + "[w, u, v] wind velocity";
      }
      return info+"[w, u] wind velocity";
    }
    return info+"vertical w wind velocity";
  }

  else // means SDnspace > 3
  {
    const std::string errmsg("SDnspace > 3 is invalid");
    throw std::invalid_argument(errmsg);
  }
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
      throw std::invalid_argument("wind velocity vectors are "
                                "not consistent with SDnspace");
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