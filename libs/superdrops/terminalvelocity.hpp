/*
 * ----- CLEO -----
 * File: terminalvelocity.hpp
 * Project: superdrops
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 8th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for terminal velocity
 * formulas (used by some types of superdroplet
 * motion and collision kernels). Formulas are 
 * contained in structures which satisfy
 * requirements of the VelocityFormula concept
*/

#ifndef TERMINALVELOCITY_HPP
#define TERMINALVELOCITY_HPP

#include <concepts>

#include "./superdrop.hpp"

template <typename V>
concept VelocityFormula = requires(V v, const Superdrop &drop)
/* Objects that are of type 'VelocityFormula'
take a superdroplet and returns something convertible
to a double (hopefully a velocity!) */
{
  {
    v(drop)
    } -> std::convertible_to<double>;
};

struct NullTerminalVelocity
{
  double operator()(const Superdrop &drop) const
  /* returns 0.0 as terminal velocity of a superdroplet */
  {
    return 0.0;
  }
};

#endif // TERMINALVELOCITY_HPP
