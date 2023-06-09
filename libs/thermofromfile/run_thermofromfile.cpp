// Author: Clara Bayley
// File: run_thermofromfile.hpp
/* Header file for functions specifically
to run uncoupled SDM where thermodynamics
are read from file */

#include "run_thermofromfile.hpp"

void receive_thermodynamics_from_thermodyn(const size_t ngbxs, 
                            const ThermodynamicsFromFile &thermodyn,
                            Kokkos::View<GridBox*> h_gridboxes)
/* Sets current thermodynamic state of SDM (time, p, temp, qv, etc.)
to match that given by the ThermodnamicsFromFile 'thermodyn' */
{
  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    auto &gbx = h_gridboxes(ii);

    gbx.state.press = thermodyn.get_press(gbx.gbxindex);
    gbx.state.temp = thermodyn.get_temp(gbx.gbxindex);
    gbx.state.qvap = thermodyn.get_qvap(gbx.gbxindex);
    gbx.state.qcond = thermodyn.get_qcond(gbx.gbxindex);

    gbx.state.wvel = thermodyn.get_wvelzfaces(gbx.gbxindex);
    gbx.state.uvel = thermodyn.get_uvelxfaces(gbx.gbxindex);
    gbx.state.vvel = thermodyn.get_vvelyfaces(gbx.gbxindex);
  }
}

void receive_thermodynamics(const int t_mdl, const int couplstep,
                            const size_t ngbxs,
                            const ThermodynamicsFromFile &thermodyn,
                            Kokkos::View<GridBox *> h_gridboxes)
/* updates time in each gbx thermodynamic state to
match t_mdl and receives thermodynamics from thermodyanmic
solver 'thermodyn' if on couplstep */
{
  const double time = step2dimlesstime(t_mdl);
  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    h_gridboxes(ii).state.time = time;
  }

  if (t_mdl % couplstep == 0)
  {
    receive_thermodynamics_from_thermodyn(ngbxs, thermodyn, h_gridboxes);
  }
}