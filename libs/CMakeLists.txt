# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# e.g. run cmake with these compiler options:
# CXX=g++ CC=gcc cmake -S ./ -B ./build

cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) # if using Kokkos with NVC++ compiler 

# essential CLEO SDM libaries #
add_subdirectory(initialisation)
add_subdirectory(superdrop_solver)
add_subdirectory(sdmgridboxes)
add_subdirectory(observers)
add_subdirectory(zarrstorage)

# libraries particular to certain CLEO configurations #
add_subdirectory(cvodecoupld EXCLUDE_FROM_ALL)
add_subdirectory(thermofromfile EXCLUDE_FROM_ALL)