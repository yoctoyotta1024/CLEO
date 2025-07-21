/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: copyfiles2txt.cpp
 * Project: configuration
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality to open files given their filenames and
 * copy their contents line by line into a .txt file.
 * Useful for copying the details of a model setup
 * e.g. configuration files and values of constants
 */

#include "configuration/copyfiles2txt.hpp"

/* open a file called filename and copy
text line by line into wfile */
void copyfile(std::ofstream &wfile, const std::filesystem::path filename);

/* creates new empty file called setup_filename and copies contents of
files listed in files2copy vector one by one */
void copyfiles2txt(const std::filesystem::path setup_filename,
                   const std::vector<std::filesystem::path> &files2copy) {
  const auto setup_filestr = setup_filename.string();
  std::cout << "----- writing to new setup file: " << setup_filestr << " -----\n";

  // Create parent directory(s) for setup_filename if not existing
  const std::filesystem::path parent_dir = setup_filename.parent_path();
  if (!std::filesystem::exists(parent_dir)) {
    std::cout << "creating directories: " << parent_dir << "\n";
    std::filesystem::create_directories(parent_dir);
  }

  std::ofstream wfile;

  wfile.open(setup_filestr, std::ios::out | std::ios::trunc);  // clear previous contents
  wfile.close();

  wfile.open(setup_filestr, std::ios::app);  // copy files one by one
  for (auto &filename : files2copy) {
    copyfile(wfile, filename);
  }
  wfile.close();

  std::cout << "---- copy complete, setup file closed -----\n";
}

/* open a file called filename and copy
text line by line into wfile */
void copyfile(std::ofstream &wfile, const std::filesystem::path filename) {
  const auto filestr = filename.string();
  std::ifstream readfile(filestr);

  std::cout << " copying " + filestr + " to setup file\n";

  wfile << "// ----------------------------- //\n";
  wfile << "// --------- " + filestr + " --------- //\n";
  wfile << "// ----------------------------- //\n";

  std::string line;
  // read file line by line
  while (getline(readfile, line)) {
    wfile << line << '\n';  // output lines to .txt file on disk
  }

  wfile << "// ----------------------------- //\n\n\n\n";
  readfile.close();
}
