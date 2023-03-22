# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build
# or CXX=gcc cmake -S ./ -B ./build
cmake_minimum_required(VERSION 3.18.0)

add_subdirectory(initialisation)
add_subdirectory(thermo_solver)
add_subdirectory(superdrop_solver)