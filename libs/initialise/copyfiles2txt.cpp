/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: copyfiles2txt.cpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
 * Modified By: CB
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

#include "initialise/copyfiles2txt.hpp"

/* open a file called filename and copy
text line by line into wfile */
void copyfile(std::ofstream &wfile, const std::string filename);

/* creates new empty file called setuptxt and copies contents of
files listed in files2copy vector one by one */
void copyfiles2txt(const std::string setuptxt, const std::vector<std::string> files2copy) {
  std::cout << "----- writing to new setuptxt file: " << setuptxt << " -----\n";

  std::ofstream wfile;

  wfile.open(setuptxt, std::ios::out | std::ios::trunc);  // clear previous contents
  wfile.close();

  wfile.open(setuptxt, std::ios::app);  // copy files one by one
  for (auto &filename : files2copy) {
    copyfile(wfile, filename);
  }
  wfile.close();

  std::cout << "---- copy complete, setuptxt file closed -----\n";
}

/* open a file called filename and copy
text line by line into wfile */
void copyfile(std::ofstream &wfile, const std::string filename) {
  std::ifstream readfile(filename);

  std::cout << " copying " + filename + " to setuptxt file\n";

  wfile << "// ----------------------------- //\n";
  wfile << "// --------- " + filename + " --------- //\n";
  wfile << "// ----------------------------- //\n";

  std::string line;
  // read file line by line
  while (getline(readfile, line)) {
    wfile << line << '\n';  // output lines to .txt file on disk
  }

  wfile << "// ----------------------------- //\n\n\n\n";
  readfile.close();
}
