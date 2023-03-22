# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build
# or CXX=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# setlibrary name
set(LIBNAME "runmodel")

# print where exactly project source dir is for this CMakeLists.txt
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

# add executables to library
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
add_library("${LIBNAME}" SHARED ${SOURCES})

# add directories to include for target
target_include_directories(${LIBNAME} PRIVATE ${PROJECT_SOURCE_DIR})

# link libraries in libs to this target
target_link_libraries("${LIBNAME}" PRIVATE thermo_solver superdrop_solver initialisation)

# some compilation properties
set_target_properties(${LIBNAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

target_compile_options(${LIBNAME} PRIVATE -Wall -g)
#target_compile_options(${LIBNAME} PRIVATE -Wall -Werror -pedantic -g)