/*
 * ----- CLEO -----
 * File: movesupersindomain.hpp
 * Project: sdmdomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Functionality related to moving superdroplets
 * (both updating their spatial coordinates and
 * moving them between gridboxes)
 */


#ifndef MOVESUPERSINDOMAIN_HPP  
#define MOVESUPERSINDOMAIN_HPP  

#include "initialise/config.hpp"
#include "initialise/timesteps.hpp"

struct MoveSupersInDomain
{
  MoveSupersInDomain(const Config &config, const Timesteps &tsteps){}
};

#endif // MOVESUPERSINDOMAIN_HPP