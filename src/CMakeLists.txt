# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# and then run cmake with this compiler:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# set project name and source directory of .cpp files
project("coupledCVODECLEO")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

add_subdirectory(superdrop_solver)
add_subdirectory(thermo_solver)

# add executables to project (only one main() function!)
set(SOURCES "main.cpp"
"run_coupledmodel.cpp"
  "./runmodel/initialisation/readbinary.cpp"
  "./runmodel/initialisation/read_gbxboundaries.cpp"
  "./runmodel/initialisation/read_initsuperdrops.cpp"
  "./runmodel/config.cpp"
  "./runmodel/gridbox.cpp"
  "./runmodel/maps4gridboxes.cpp"
  "./runmodel/observers/observers.cpp"
  "./runmodel/observers/observer_singlevariable.cpp"
  "./runmodel/observers/observer_thermostate.cpp"
  "./runmodel/observers/zarrstores.cpp"
  "./runmodel/superdrops_with_gridboxes.cpp"
  "./runmodel/timesteps.cpp"
  )

add_executable("${PROJECT_NAME}" ${SOURCES})

target_link_libraries("${PROJECT_NAME}" PUBLIC thermo_solver)
target_link_libraries("${PROJECT_NAME}" PUBLIC superdrop_solver)

# # add directories to include for targets to us
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS OFF
  CXX_STANDARD 20)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -g)
#target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -pedantic -g)