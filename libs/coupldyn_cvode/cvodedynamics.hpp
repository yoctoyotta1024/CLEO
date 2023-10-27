/*
 * ----- CLEO -----
 * File: cvodedynamics.hpp
 * Project: coupldyn_cvode
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 27th October 2023
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
 * two-way to cvode adiabatic parcel ODE solver
 */

#ifndef CVODEDYNAMICS_HPP 
#define CVODEDYNAMICS_HPP 

#include <iostream>
#include <vector>

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */

#include "../cleoconstants.hpp"
#include "./differential_functions.hpp"
#include "initialise/config.hpp"

namespace dlc = dimless_constants;

struct CvodeDynamics
/* type satisfying CoupledDyanmics solver
concept specifically for thermodynamics
of adiabatically expanding parcel (0-D) */
{
private:
  const unsigned int interval;
  std::vector<double> previousstates; // holds states press, temp, qvap and qcond before timestep iterated

  /* SUNDIALS CVODE solver stuff */
  SUNContext sunctx;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;

  /* ODE problem stuff */
  static constexpr int NVARS = 4;    // no. of distinct variables (= no. ODEs per grid box)
  const size_t neq;              // No. of equations/ODEs (= total no. variables across all Grid Boxes)
  realtype t;
  realtype RTOL;
  N_Vector y;
  N_Vector re_y;
  N_Vector ATOLS;
  UserData data;

  int check_retval(void *returnvalue, const char *funcname, int opt);
  /* Check function return value for memory or sundials CVODE error */

  void run_dynamics(const unsigned int t_mdl,
                    const unsigned int t_next) const;
  
public:
  CvodeDynamics(const Config &config,
                const unsigned int couplstep)
      : interval(couplstep)
  {
    /* CVODE thermodynamics solver */
    CvodeThermoSolver cvode(config,
                            initcvodethermo(sdm.gbxmaps.ngridboxes, config));
  }

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
      run_dynamics(t_mdl, t_next);
    }
  }

};

#endif // CVODEDYNAMICS_HPP