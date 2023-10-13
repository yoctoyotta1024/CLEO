/*
 * ----- CLEO -----
 * File: gridbox.hpp
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
 * Functions and structures related to the CLEO gridboxes
 */


#ifndef GRIDBOX_HPP 
#define GRIDBOX_HPP 

struct Gridbox
/* gridbox contains vector of superdroplets in grid box,
thermodynamic state temp, pressure, etc. used for SDM,
and index for finding associated grridbox in
coupled thermodynamics */
{
  unsigned int gbxindex; // index (unique identifier) of gridbox
  std::shared_ptr<Detectors> detectors;
  std::span<SuperdropWithGbxindex> span4SDsinGBx;
  ThermoState state;
};

struct Gridboxes
{
  // kokkos array (dualview)
};

#endif // GRIDBOX_HPP