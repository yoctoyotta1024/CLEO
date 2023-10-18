/*
 * ----- CLEO -----
 * File: sortsupers.hpp
 * Project: gridboxes
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions used when sorting /shuffling superdrops
 * e.g. based on their gridbox indexes
 */


#ifndef SORTSUPERS_HPP
#define SORTSUPERS_HPP

#include <iostream>

#include "../kokkosaliases.hpp"
#include "superdrops/superdrop.hpp"

viewd_supers sort_supers(viewd_supers supers);
/* sort a view of superdroplets by their sdgbxindexes */

#endif // SORTSUPERS_HPP