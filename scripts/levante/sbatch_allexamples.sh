#!/bin/bash

path2CLEO=${1:-${HOME}/CLEO}

if [ -z "${path2CLEO}" ];
then
  echo "please specify the path2CLEO"
else
  sbatch ${path2CLEO}/examples/adiabaticparcel/as2017.sh
  sbatch ${path2CLEO}/examples/adiabaticparcel/cuspbifurc.sh
  sbatch ${path2CLEO}/examples/boxmodelcollisions/shima2009.sh
  sbatch ${path2CLEO}/examples/boxmodelcollisions/breakup.sh
  #sbatch ${path2CLEO}/examples/bubble3d/bubble3d.sh # TODO(CB): complete
  sbatch ${path2CLEO}/examples/constthermo2d/constthermo2d.sh
  sbatch ${path2CLEO}/examples/divfreemotion/divfree2d.sh
  sbatch ${path2CLEO}/examples/eurec4a1d/eurec4a1d.sh
  sbatch ${path2CLEO}/examples/fromfile/fromfile.sh
  sbatch ${path2CLEO}/examples/fromfile_irreg/fromfile_irreg.sh
  sbatch ${path2CLEO}/examples/rainshaft1d/rainshaft1d.sh
  sbatch ${path2CLEO}/examples/speedtest/speedtest.sh
fi
