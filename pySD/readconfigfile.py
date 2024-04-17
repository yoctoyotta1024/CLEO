'''
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: readconfigfile.py
Project: pySD
Created Date: Wednesday 17th April 2024
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

def read_configparams_into_floats(filename):
  """returns dictionary of value: float from
  values assigned in a config .txt file.
  Also returns dictionary of notfloats
  for values that couldn't be converted. """

  floats = {}
  notfloats = {}
  with open(filename) as file:
    rlines=[]
    filelines = file.readlines()
    for line in filelines:
      if line[0] != "#" and line[0] != "/" and "=" in line:
        goodline = remove_excess_line(line)
        rlines.append(goodline)

    for line in rlines:
      ind = line.find("=")
      name =  line[:ind]
      value =  line[ind+1:]

      try:
        floats[name] = float(value)
      except ValueError:
        notfloats[name] = value

  try:
    floats["nspacedims"] = int(floats["nspacedims"])              # no spatial coords to SDs
  except:
    pass

  return floats
