# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build
# or CXX=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# set project name and source directory of .cpp files
project("testing")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# add executables to project (only one main() function!)
set(SOURCES "testing.cpp")
add_executable("${PROJECT_NAME}" ${SOURCES})

# # add directories to include for targets to us
# set(CLEODIR "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO")
set(CLEODIR "../")
include_directories("${CLEODIR}/src/runmodel")
include_directories("${CLEODIR}/libs")
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CLEODIR})

set(LIBNAME "runmodel")
# add executables to runmodel library
set(LIBSOURCES 
    "../src/runmodel/gridbox.cpp"
    "../src/runmodel/maps4gridboxes.cpp"
    "../src/runmodel/observers/observers.cpp"
    "../src/runmodel/observers/observer_moments.cpp"
    "../src/runmodel/observers/observer_singlevariable.cpp"
    "../src/runmodel/observers/observer_thermostate.cpp"
    "../src/runmodel/observers/zarrstores.cpp"
    "../src/runmodel/superdrops_in_gridboxes.cpp"
    "../src/runmodel/timesteps.cpp"
  )
add_library("${LIBNAME}" SHARED ${LIBSOURCES})

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/initialisation" ${CMAKE_CURRENT_BINARY_DIR}/libs/initialisation)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/superdrop_solver" ${CMAKE_CURRENT_BINARY_DIR}/libs/superdrop_solver)

target_link_libraries("${LIBNAME}" PRIVATE superdrop_solver initialisation)
target_link_libraries("${PROJECT_NAME}" PUBLIC superdrop_solver initialisation runmodel)

set_target_properties(${LIBNAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

set_target_properties(${PROJECT_NAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

target_compile_options(${LIBNAME} PUBLIC -Wall -pedantic -g -gdwarf-4)
target_compile_options(${PROJECT_NAME} PUBLIC -Wall -pedantic -g -gdwarf-4)
#target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -pedantic -g)