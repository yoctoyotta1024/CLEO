/*
 * ----- CLEO -----
 * File: fromfiledynamics.hpp
 * Project: coupldyn_fromfile
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 6th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct obeying coupleddyanmics concept for
 * dynamics solver in CLEO where coupling is
 * one-way and dynamics are read from file
 */

#ifndef FROMFILEDYNAMICS_HPP 
#define FROMFILEDYNAMICS_HPP 

#include <array>
#include <vector>
#include <memory>

#include "initialise/config.hpp"

struct DynamicsVariables
/* contains 1-D vector for each (thermo)dynamic
variable which is ordered by gridbox at every timestep
e.g. press = [p_gbx0(t0), p_gbx1(t0), ,... , p_gbxN(t0), 
p_gbx0(t1), p_gbx1(t1), ..., p_gbxN(t1), ..., p_gbxN(t_end)]
"pos[_X]" gives position of variable in a vector to read
current timestep from for the first gridbox (gbx0)  */
{
  /* (thermo)dynamic variables read from file */
  std::vector<double> press;
  std::vector<double> temp;
  std::vector<double> qvap;
  std::vector<double> qcond;
  std::vector<double> wvel; // w velocity define of z faces of gridboxes
  std::vector<double> uvel; // u velocity define of x faces of gridboxes
  std::vector<double> vvel; // v velocity define of y faces of gridboxes

  /* position in vector for 0th gridbox at current timestep  */
  const std::array<size_t, 3> ndims; // number of (centres of) gridboxes in [z,x,y] directions
  size_t pos;                        // for variable defined at gridbox centres
  size_t pos_zface;                  // for variable defined at gridbox z faces
  size_t pos_xface;                  // for variable defined at gridbox x faces
  size_t pos_yface;                  // for variable defined at gridbox y faces

  DynamicsVariables(const Config &config);

  void increment_position();
  /* updates positions to gbx0 in vector (for
  acessing value at next timestep) */
};


struct FromFileDynamics
/* type satisfying CoupledDyanmics solver concept
specifically for thermodynamics and wind velocities
that are read from binary files */
{
private:
  const unsigned int interval;
  std::unique_ptr<DynamicsVariables> dynvars; // pointer to (thermo)dyanmic variables

  void run_dynamics(const unsigned int t_mdl) const
  /* increment position of thermodata for 0th gridbox
  to positon at next timestep (ie. ngridbox_faces
  further along vector) */
  {
    dynvars->increment_position();
  }

public:
  FromFileDynamics(const Config &config,
                   const unsigned int couplstep)
      : interval(couplstep),
        dynvars(std::make_unique<DynamicsVariables>(config)) {}

  auto get_couplstep() const
  {
    return interval;
  }

  void prepare_to_timestep() const;

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void run_step(const unsigned int t_mdl,
                const unsigned int t_next) const
  {
    if (on_step(t_mdl))
    {
      run_dynamics(t_mdl);
    }
  }

};

#endif // FROMFILEDYNAMICS_HPP