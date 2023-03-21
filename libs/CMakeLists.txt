# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build
# or CXX=gcc cmake -S ./ -B ./build
cmake_minimum_required(VERSION 3.18.0)

project(project
  LANGUAGES CXX
  DESCRIPTION "CLEO superdroplet and CVODE libraries by Clara Bayley"
  VERSION 1.0
)

add_subdirectory(thermo_solver)
add_subdirectory(superdrop_solver)
add_subdirectory(initialisation)