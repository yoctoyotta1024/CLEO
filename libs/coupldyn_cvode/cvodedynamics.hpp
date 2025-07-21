/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cvodedynamics.hpp
 * Project: coupldyn_cvode
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct obeying coupleddynamics concept for dynamics solver in CLEO where coupling is two-way
 * to cvode adiabatic parcel ODE solver
 */

#ifndef LIBS_COUPLDYN_CVODE_CVODEDYNAMICS_HPP_
#define LIBS_COUPLDYN_CVODE_CVODEDYNAMICS_HPP_

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../cleoconstants.hpp"
#include "configuration/optional_config_params.hpp"
#include "coupldyn_cvode/differentialfuncs.hpp"

namespace dlc = dimless_constants;

/* type satisfying CoupledDyanmics solver
concept specifically for thermodynamics
of adiabatically expanding parcel (0-D) */
struct CvodeDynamics {
 private:
  const unsigned int interval;
  std::function<double(unsigned int)>
      step2dimlesstime;  // function to convert timesteps to real time

  /* SUNDIALS CVODE solver stuff */
  SUNContext sunctx;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;

  /* ODE problem stuff */
  static constexpr int NVARS = 4;  // no. of distinct variables (= no. ODEs per grid box)
  const size_t neq;  // No. of equations/ODEs (= total no. variables across all Grid Boxes)
  realtype t;
  realtype RTOL;
  N_Vector y;
  N_Vector re_y;
  N_Vector ATOLS;
  UserData data;
  std::vector<double>
      previousstates;  // holds states press, temp, qvap and qcond before timestep iterated

  /* print initial ODE setup to the terminal screen */
  void print_initODEstatement() const;

  int run_dynamics(const unsigned int t_next);

  /* return vector of dimensionless initial conditions
  for thermodynamic variables (p, temp, qv, qc) to
  initialise cvode thermodynamics solver */
  std::vector<double> initial_conditions(
      const OptionalConfigParams::CvodeDynamicsParams &config) const;

  /* set values in UserData structure for odes_func */
  void init_userdata(const size_t neq, const double wmax, const double tauhalf);

  /* function does all the setup steps in order
  to use CVODE sundials ODE solver */
  int setup_ODE_solver(const double i_rtol, const double i_atol);

  /* Check function return value for memory or sundials CVODE error */
  int check_retval(void *returnvalue, const char *funcname, int opt);

 public:
  /* construct instance of CVODE ODE
  solver with initial conditions */
  CvodeDynamics(const OptionalConfigParams::CvodeDynamicsParams &config,
                const unsigned int couplstep,
                const std::function<double(unsigned int)> step2dimlesstime);

  /* print final statistics to the
  terminal screen and free CVODE memory */
  ~CvodeDynamics();

  auto get_couplstep() const { return interval; }

  double get_time() const { return t; }

  /* returns ii'th previous state [press, temp, qvap, qcond] */
  std::array<double, 4> get_previous_state(const size_t ii) const {
    const size_t jj(NVARS * ii);

    std::array<double, 4> prevstate;
    // n = 0,1,2,3
    for (size_t n(0); n < 4; ++n) {
      prevstate.at(n) = previousstates.at(jj + n);
    }

    return prevstate;
  }

  /* returns ii'th [press, temp, qvap, qcond] state */
  std::array<double, 4> get_current_state(const size_t ii) const {
    const size_t jj(NVARS * ii);

    std::array<double, 4> state;
    // n = 0,1,2,3
    for (size_t n(0); n < 4; ++n) {
      state.at(n) = NV_Ith_S(y, jj + n);  // state
    }

    return state;
  }

  /* Reinitialize the solver after discontinuous change
  in temp, qv and qc (e.g. due to condensation) */
  int reinitialise(const double next_t, const std::vector<double> &delta_y);

  /* checks initial y has been set and then
  prints statement about cvode ODEs configuration */
  void prepare_to_timestep() const;

  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  void run_step(const unsigned int t_mdl, const unsigned int t_next) {
    if (on_step(t_mdl)) {
      run_dynamics(t_next);
    }
  }
};

#endif  // LIBS_COUPLDYN_CVODE_CVODEDYNAMICS_HPP_
