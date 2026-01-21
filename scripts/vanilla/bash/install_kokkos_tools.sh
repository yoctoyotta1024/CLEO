#!/bin/bash

### ------------------------------------------------------- ###
### running script sucessfully installs kokkostools for a gcc compiler with
### openmpi on a "vanilla" computer assuming you have already performed
### ```git clone git@github.com:kokkos/kokkos-tools.git```
### in ${path2toolsrepo}
### ------------------------------------------------------- ###

path2toolsrepo=$1 # absolute path for kokkos-tools repository clone installation
compilername=${2:-gcc} # compile tools with "gcc"
root4tools=${3:-${CLEO_KOKKOSTOOLS}} # absolute path for kokkos-tools installation

if [ "${compilername}" == "" ]
then
  echo "Bad input, please specify compiler name to build kokkos-tools with"
  exit 1
elif [[ "${compilername}" != "gcc" ]]
then
  echo "Bad input, unrecognised compiler name"
  exit 1
fi

if [[ "${path2toolsrepo}" == "" || "${root4tools}" == "" ]]
then
  echo "Bad inputs, please specify absolute path for kokkos-tools repo and where you want to install the tools"
  exit 1
else
  CXX="$(command -v mpic++)"

  cd ${path2toolsrepo}
  if [ ! -d "./kokkos-tools" ]
  then
    echo "ERROR: kokkos-tools source directory not found. Please clone it into path2toolsrepo" >&2
    exit 1
  fi

  mkdir ${root4tools}
  cd ./kokkos-tools && mkdir ./myBuild
  cmake -DCMAKE_CXX_COMPILER=${CXX} \
      -S ./ -B ./myBuild \
      -DCMAKE_INSTALL_PREFIX=${root4tools}
  cd ./myBuild && make && make install

  echo "----- INSTALLATION SUMMARY -----"
  echo "in ${root4tools}/bin:" && ls ${root4tools}/bin
  if [ -d "${root4tools}/lib" ]
  then
    echo "in ${root4tools}/lib:" && ls ${root4tools}/lib
  fi
  if [ -d "${root4tools}/lib64" ]
  then
    echo "in ${root4tools}/lib64:" && ls ${root4tools}/lib64
  fi
  echo "--------------------------------"
  echo "SUCCESS: kokkos-tools installed in ${root4tools}"
fi

### ------------ Notes on using profiler for executable ------------ ###
# example for for tools installed in /path/to/tools/kokkostools/ on macOS:
# A) see tool libraries installed in /path/to/tools/kokkostools/lib/
# B) export required tool library, e.g.
#     e.g. export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib/libkp_kernel_timer.dylib
#      or  export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib/libkp_space_time_stack.dylib
# C) run executable ./[exec].exe (kokkos initialise loads dynamic library pointers)
# D) read *.dat output
#     e.g. with kp reader
#          export DYLD_LIBRARY_PATH=/path/to/tools/kokkostools/lib/:$LD_LIBRARY_PATH
#          /path/to/tools/kokkostools/bin/kp_reader *.dat > ./bin/kp_kernel_timer.txt
#    or pipe kp_space_time_stack output durign runtime: ./[exec].exe > runtime_output.txt
# E) Also note useful debugging tool to find where program crashed (e.g. inside kernel):
#   export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib/libkp_kernel_logger.dylib
### ---------------------------------------------------------------- ###
