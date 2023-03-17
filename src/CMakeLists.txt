# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# or enter build dir and just do " cmake [source dir] "

#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.23.1")

# add compiler flags (TOD0: do this in proper way...)
add_compile_options(-Wall -g -std=c++20)

#add_compile_options(-Wall -Werror -pedantic -g -std=c++20)

# # specify the C++ standard (or use flag in line 8 above)
# set(CMAKE_CXX_STANDARD 11)
# target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)
# set_target_properties(${PROJECT_NAME} PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON)
# set_target_properties(${PROJECT_NAME} PROPERTIES CMAKE_CXX_EXTENSIONS ON)

# set project name and source directory of .cpp files
project("coupledmodel")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# add executables to project (only one main() function!)
set(SOURCES "main.cpp"
  "./thermodynamic_equations.cpp"
  "./coupledmodel/initialisation/readbinary.cpp"
  "./coupledmodel/initialisation/read_gbxboundaries.cpp"
  "./coupledmodel/initialisation/read_initsuperdrops.cpp"
  "./coupledmodel/config.cpp"
  "./coupledmodel/gridbox.cpp"
  "./coupledmodel/maps4gridboxes.cpp"
  "./coupledmodel/observers/observers.cpp"
  "./coupledmodel/observers/observer_singlevariable.cpp"
  "./coupledmodel/observers/observer_thermostate.cpp"
  "./coupledmodel/observers/zarrstores.cpp"
  "./coupledmodel/run_coupledmodel.cpp"
  "./coupledmodel/superdrops_with_gridboxes.cpp"
  "./coupledmodel/timesteps.cpp"
  "./thermo_solver/cvodethermosolver.cpp"
  "./thermo_solver/differential_functions.cpp"
  "./superdrop_solver/condensationmethod.cpp"
  "./superdrop_solver/impliciteuler.cpp"
  "./superdrop_solver/superdrop.cpp")
add_executable("${PROJECT_NAME}" ${SOURCES})

# # add directories to include for targets to us
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR})

# CVODE libaries stuff
# specify libraries to link against
set(LIBDIR /usr/local/cvode-6/lib)
find_library(CVODE_LIB NAMES sundials_cvode HINTS "${LIBDIR}")
find_library(NVECS_LIB NAMES sundials_nvecserial HINTS "${LIBDIR}")

# find_library(MANYVECS_LIB NAMES sundials_nvecmanyvector  HINTS "${LIBDIR}")
set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB})
target_link_libraries(${PROJECT_NAME} PRIVATE ${SUNDIALS_LIBS})

# CVODE include libraries stuff
target_include_directories(${PROJECT_NAME} PRIVATE /usr/local/cvode-6/include)

# ## some extra stuff that's redundant at the moment
# include_directories(${PROJECT_SOURCE_DIR})
# install(TARGETS "${PROJECT_NAME}" DESTINATION bin) # ${CMAKE_INSTALL_PREFIX}/bin
# install(FILES "testing.cpp" DESTINATION src) # ${CMAKE_INSTALL_PREFIX}/src
