/*
 * ----- CLEO -----
 * File: contigraggedstorage.hpp
 * Project: zarr
 * Created Date: Monday 23rd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 23rd October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * File for Contiguous Ragged Array Storage
 * used to store superdroplet attributes
 * (see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
 * in a FFStore obeying zarr storage specification verion 2:
 * https://zarr.readthedocs.io/en/stable/spec/v2.html */


#ifndef CONTIGRAGGEDSTORAGE
#define CONTIGRAGGEDSTORAGE 

#include <concepts>
#include <vector>
#include <string>
#include <utility>
#include <tuple>

#include "./fsstore.hpp"
#include "./storehelpers.hpp"
#include "./superdropbuffers.hpp"
#include "../cleoconstants.hpp"
#include "superdrops/superdrop.hpp"

template <SuperdropBuffers SDIntoStore>
class ContigRaggedStorage
/* Class for outputting Superdrop's data into zarr storage in
arrays of contigous ragged representation with 'chunkcount' number
of chunks that have a fixed chunksize. Works by filling buffers in
sdbuffers with superdrop data and then writing these buffers
into chunks in their corresponding array stores when number of
datapoints copied to the buffers reaches chunksize. */
{
}

#endif // CONTIGRAGGEDSTORAGE