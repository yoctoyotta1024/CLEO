/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: cvodedynamics.hpp
 * Project: coupldyn_cvode
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 21st November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct obeying coupleddynamics concept for
 * dynamics solver in CLEO where coupling is
 * two-way to cvode adiabatic parcel ODE solver
 */

#ifndef CVODEDYNAMICS_HPP
#define CVODEDYNAMICS_HPP

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */

#include "../cleoconstants.hpp"
#include "./differentialfuncs.hpp"
#include "initialise/config.hpp"

namespace dlc = dimless_constants;

struct CvodeDynamics
/* type satisfying CoupledDyanmics solver
concept specifically for thermodynamics
of adiabatically expanding parcel (0-D) */
{
private:
  const unsigned int interval;
  std::function<double(unsigned int)> step2dimlesstime; // function to convert timesteps to real time

  /* SUNDIALS CVODE solver stuff */
  SUNContext sunctx;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;

  /* ODE problem stuff */
  static constexpr int NVARS = 4;    // no. of distinct variables (= no. ODEs per grid box)
  const size_t neq;                  // No. of equations/ODEs (= total no. variables across all Grid Boxes)
  realtype t;
  realtype RTOL;
  N_Vector y;
  N_Vector re_y;
  N_Vector ATOLS;
  UserData data;
  std::vector<double> previousstates; // holds states press, temp, qvap and qcond before timestep iterated

  void print_initODEstatement() const;
  /* print initial ODE setup to the terminal screen */

  int run_dynamics(const unsigned int t_next);

  std::vector<double> initial_conditions(const Config &config) const;
  /* return vector of dimensionless initial conditions
  for thermodynamic variables (p, temp, qv, qc) to
  initialise cvode thermodynamics solver */

  void init_userdata(const size_t neq,
                     const double wmax,
                     const double tauhalf);
  /* set values in UserData structure for odes_func */

  int setup_ODE_solver(const double i_rtol, const double i_atol);
  /* function does all the setup steps in order
  to use CVODE sundials ODE solver */

  int check_retval(void *returnvalue, const char *funcname, int opt);
  /* Check function return value for memory or sundials CVODE error */

public:
  CvodeDynamics(const Config &config,
                const unsigned int couplstep,
                const std::function<double(unsigned int)> step2dimlesstime);
  /* construct instance of CVODE ODE
  solver with initial conditions */

  ~CvodeDynamics();
  /* print final statistics to the
  terminal screen and free CVODE memory */

  auto get_couplstep() const { return interval; }

  double get_time() const { return t; }

  std::array<double, 4> get_previous_state(const size_t ii) const
  /* returns ii'th previous state [press, temp, qvap, qcond] */
  {
    const size_t jj(NVARS * ii);

    std::array<double, 4> prevstate;
    for (size_t n(0); n < 4; ++n) // n = 0,1,2,3
    {
      prevstate.at(n) = previousstates.at(jj + n);
    }

    return prevstate;
  }

  std::array<double, 4> get_current_state(const size_t ii) const
  /* returns ii'th [press, temp, qvap, qcond] state */
  {
    const size_t jj(NVARS * ii);

    std::array<double, 4> state;
    for (size_t n(0); n < 4; ++n) // n = 0,1,2,3
    {
      state.at(n) = NV_Ith_S(y, jj + n); // state
    }

    return state;
  }

  int reinitialise(const double next_t,
                   const std::vector<double> &delta_y);
  /* Reinitialize the solver after discontinuous change
  in temp, qv and qc (e.g. due to condensation) */

  void prepare_to_timestep() const;
  /* checks initial y has been set and then
  prints statement about cvode ODEs configuration */

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void run_step(const unsigned int t_mdl,
                const unsigned int t_next)
  {
    if (on_step(t_mdl))
    {
      run_dynamics(t_next);
    }
  }

};

#endif // CVODEDYNAMICS_HPP
