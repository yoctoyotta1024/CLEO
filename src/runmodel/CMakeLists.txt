# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build
# or CXX=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# set project name and source directory of .cpp files
project("runmodel")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

# add executables to project (only one main() function!)
set(SOURCES 
    "gridbox.cpp"
    "maps4gridboxes.cpp"
    "observers/observers.cpp"
    "observers/observer_singlevariable.cpp"
    "observers/observer_thermostate.cpp"
    "observers/zarrstores.cpp"
    "superdrops_with_gridboxes.cpp"
    "timesteps.cpp"
  )
add_library("${PROJECT_NAME}" SHARED ${SOURCES})

# # add directories to include for targets to us
target_link_libraries("${PROJECT_NAME}" PRIVATE thermo_solver superdrop_solver initialisation)

# # add directories to include for targets to us
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -g)
#target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -pedantic -g)