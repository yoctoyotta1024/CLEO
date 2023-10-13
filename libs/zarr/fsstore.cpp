/*
 * ----- CLEO -----
 * File: fsstore.cpp
 * Project: zarr
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
 * some functionality for using stores obeying the
 * zarr storage specification version 2 (e.g. see FSStore)
 * https://zarr.readthedocs.io/en/stable/spec/v2.html
 */


#include "./fsstore.hpp"

bool FSStore::write(std::string_view key,
                    std::span<const uint8_t> buffer)
/* write function called by StoreAccessor
once data has been converted into a vector
of unsigned integer types */
{
  auto path = basedir / key;
  auto mode = std::ios::out | std::ios::binary;
  std::ofstream out(path, mode);
  if (!out.good())
  {
    std::cout << "couldn't open " << path << ",\n "
              << "making directory " << path.parent_path() << "\n";
    std::filesystem::create_directories(path.parent_path());
    out.open(path, mode);
  }
  if (!out.good())
  {
    std::cout << "can't write to " << path << "\n";
    return false;
  }
  out.write(reinterpret_cast<const char *>(buffer.data()), buffer.size());
  return true;
}