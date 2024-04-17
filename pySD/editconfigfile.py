'''
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: editconfigfile.py
Project: pySD
Created Date: Wednesday 17th January 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 17th April 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
'''

# TODO(CB): rewrite funciton for .yaml not .txt config file

def edit_config_params(filename, params2change):
  """rewrites config file with parameters listed in
  dict params2change edited to new values also in dict"""

  wlines=[]

  with open(filename) as file:
    filelines = file.readlines()
    for line in filelines:
      wlines.append(line)

  for l, line in enumerate(wlines):
    if line[0] != "#" and line[0] != "/" and "=" in line:
      for key, value in params2change.items():
        if key in line:

          # create line with new value for key
          newline = key+" = "+str(value)
          newline = newline.ljust(40)

          # add comment to new line if there is one
          ind = line.find("#")
          newline = newline+line[ind:]

          # overwrite line with newline
          wlines[l] = newline

  file = open(filename, "w")
  file.writelines(wlines)
  file.close()
