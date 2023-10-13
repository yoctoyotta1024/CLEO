/*
 * ----- CLEO -----
 * File: copyfiles2txt.hpp
 * Project: initialise
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
 *  structure for opening files given their filenames and copying their
 * contents line by line into a .txt file. Useful for copying
 * the details of a model setup e.g. configuration and
 * values of constants
 */


#ifndef COPYFILES2TXT_HPP
#define COPYFILES2TXT_HPP

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <tuple>
#include <sstream>

struct CopyFiles2Txt
{
  void operator()(const std::string setup_txt,
                  const std::vector<std::string> files2copy) const
  /* creates new empty file called setup_txt and copies contents of
  files listed in files2copy vector one by one */
  {
    std::ofstream txtfile;

    txtfile.open(setup_txt, std::ios::out | std::ios::trunc); // clear previous contents
    txtfile.close();

    txtfile.open(setup_txt, std::ios::app); // copy files one by one
    for (auto &filename : files2copy)
    {
      writefile2txt(txtfile, filename);
    }
    txtfile.close();
  }

  void writefile2txt(std::ofstream &wfile,
                     const std::string filename) const
  /* open a file called filename and copy text
  line by line into wfile */
  {
    std::ifstream readfile(filename);

    std::cout << " writing " + filename + " to setup_txt file\n";

    wfile << "// ----------------------------- //\n";
    wfile << "// --------- " + filename + " --------- //\n";
    wfile << "// ----------------------------- //\n";

    std::string line;
    while (getline(readfile, line)) // read file line by line
    {
      wfile << line << '\n'; // output lines to .txt file on disk
    }

    wfile << "// ----------------------------- //\n\n\n\n";
    readfile.close();
  }
};

#endif // COPYFILES2TXT_HPP