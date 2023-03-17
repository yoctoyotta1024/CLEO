# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# or enter build dir and just do " cmake [source dir] "

#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.23.1")

# add compiler flags (TOD0: do this in proper way...)
add_compile_options(-Wall -g -std=c++20)

# set project name and source directory of .cpp files
project("testing")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# add executables to project (only one main() function!)
set(SOURCES "testing.cpp"
    "../src/thermodynamic_equations.cpp"
    "../src/coupledmodel/initialisation/readbinary.cpp"
    "../src/coupledmodel/initialisation/read_gbxboundaries.cpp"
    "../src/coupledmodel/initialisation/read_initsuperdrops.cpp"
    "../src/coupledmodel/config.cpp"
    "../src/coupledmodel/gridbox.cpp"
    "../src/coupledmodel/maps4gridboxes.cpp"
    "../src/coupledmodel/observers/observers.cpp"
    "../src/coupledmodel/observers/observer_singlevariable.cpp"
    "../src/coupledmodel/observers/observer_thermostate.cpp"
    "../src/coupledmodel/observers/zarrstores.cpp"
    #"../src/coupledmodel/run_coupledmodel.cpp"
    "../src/coupledmodel/superdrops_with_gridboxes.cpp"
    "../src/coupledmodel/timesteps.cpp"
    #"../src/thermo_solver/cvodethermosolver.cpp"
    #"../src/thermo_solver/differential_functions.cpp"
    "../src/superdrop_solver/condensationmethod.cpp"
    "../src/superdrop_solver/impliciteuler.cpp"
    "../src/superdrop_solver/superdrop.cpp")
add_executable("${PROJECT_NAME}" ${SOURCES})


# # add directories to include for targets to us
set(INCLUDEDIR /Users/yoctoyotta1024/Documents/autumnwinter2022_23/clara-bayley-superdroplets/superdroplet_model/src)
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDEDIR})

# # CVODE libaries stuff
# # specify libraries to link against
# set(LIBDIR /usr/local/cvode-6/lib)
# find_library(CVODE_LIB NAMES sundials_cvode HINTS "${LIBDIR}")
# find_library(NVECS_LIB NAMES sundials_nvecserial HINTS "${LIBDIR}")

# # find_library(MANYVECS_LIB NAMES sundials_nvecmanyvector  HINTS "${LIBDIR}")
# set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB})
# target_link_libraries(${PROJECT_NAME} PRIVATE ${SUNDIALS_LIBS})

# # CVODE include libraries stuff
# target_include_directories(${PROJECT_NAME} PRIVATE /usr/local/cvode-6/include)
