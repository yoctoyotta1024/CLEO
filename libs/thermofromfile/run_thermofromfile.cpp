// Author: Clara Bayley
// File: run_thermofromfile.hpp
/* Header file for functions specifically
to run uncoupled SDM where thermodynamics
are read from file */

#include "run_thermofromfile.hpp"

void recieve_windcomponents(const ThermodynamicsFromFile &thermodyn,
                            std::vector<GridBox> &gridboxes)
{
  for (auto &gbx : gridboxes)
  {
    if (SDnspace >= 1)
    {
      gbx.state.wvel = thermodyn.get_wvel(gbx.gbxindex);
    }
    if (SDnspace >= 2)
    {
      gbx.state.uvel = thermodyn.get_uvel(gbx.gbxindex);
    }
    if (SDnspace >=3)
    {
      gbx.state.vvel = thermodyn.get_vvel(gbx.gbxindex);
    }
  }
}

void recieve_thermodynamics(const int time,
                            const ThermodynamicsFromFile &thermodyn,
                            std::vector<GridBox> &gridboxes)
/* Sets current thermodynamic state of SDM (time, p, temp, qv, etc.)
to match that given by the ThermodnamicsFromFile 'thermodyn' */
{
  for (auto &gbx : gridboxes)
  {
    gbx.state.time = time; 
    
    gbx.state.press = thermodyn.get_press(gbx.gbxindex);
    gbx.state.temp = thermodyn.get_temp(gbx.gbxindex);
    gbx.state.qvap = thermodyn.get_qvap(gbx.gbxindex);
    gbx.state.qcond = thermodyn.get_qcond(gbx.gbxindex);
  }  
  
  recieve_windcomponents(thermodyn, gridboxes);
}