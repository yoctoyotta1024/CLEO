/*
 * ----- CLEO -----
 * File: cartesianmotion.hpp
 * Project: cartesiandomain
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
 * Motion of a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * the sdgbxindex updated accordingly for a
 * cartesian domain with finite/periodi boundary 
 * conditions
 */


#ifndef CARTESIANMOTION_HPP
#define CARTESIANMOTION_HPP

// inline Motion<CartesianMaps> auto
// create_motion(const unsigned int motionstep)
// {
//   // using TerminalVelocity = NullTerminalVelocity;
//   // using TerminalVelocity = RogersYauTerminalVelocity;
//   // using TerminalVelocity = SimmelTerminalVelocity;

//   // return PredCorrMotion<CartesianMaps,
//   //                       TerminalVelocity>(motionstep,
//   //                                         &step2dimlesstime);
//   return NullMotion{};                                                                               
// }

#endif // CARTESIANMOTION_HPP