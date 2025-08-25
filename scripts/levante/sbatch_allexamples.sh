#!/bin/bash

path2CLEO=${1:-${HOME}/CLEO}
path2examplesbash=${path2CLEO}/scripts/levante/examples

if [ -z "${path2CLEO}" ];
then
  echo "please specify the path2CLEO"
else
  sbatch ${path2examplesbash}/as2017.sh
  sbatch ${path2examplesbash}/cuspbifurc.sh
  sbatch ${path2examplesbash}/shima2009.sh
  sbatch ${path2examplesbash}/breakup.sh
  sbatch ${path2examplesbash}/bubble3d.sh
  sbatch ${path2examplesbash}/constthermo2d.sh
  sbatch ${path2examplesbash}/divfree2d.sh
  sbatch ${path2examplesbash}/eurec4a1d.sh
  sbatch ${path2examplesbash}/fromfile.sh
  sbatch ${path2examplesbash}/fromfile_irreg.sh
  sbatch ${path2examplesbash}/python_bindings.sh
  sbatch ${path2examplesbash}/rainshaft1d.sh
  sbatch ${path2examplesbash}/kokkostools.sh
fi
