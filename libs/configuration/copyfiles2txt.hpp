/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: copyfiles2txt.hpp
 * Project: configuration
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * function(s) to open files given their filenames and
 * copy their contents line by line into a .txt file.
 * Useful for copying the details of a model setup
 * e.g. configuration files and values of constants
 */

#ifndef LIBS_CONFIGURATION_COPYFILES2TXT_HPP_
#define LIBS_CONFIGURATION_COPYFILES2TXT_HPP_

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/* creates new empty file called setup_filename and copies contents
of files listed in files2copy vector one by one */
void copyfiles2txt(const std::filesystem::path setup_filename,
                   const std::vector<std::filesystem::path> &files2copy);

#endif  // LIBS_CONFIGURATION_COPYFILES2TXT_HPP_
