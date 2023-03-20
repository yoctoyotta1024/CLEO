# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# and then run cmake with this compiler:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# set project name and source directory of .cpp files
project("coupledCVODECLEO")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# add executables to project (only one main() function!)
add_subdirectory(runmodel)
set(SOURCES "main.cpp" "run_coupledmodel.cpp")
add_executable("${PROJECT_NAME}" ${SOURCES})

set(LIBSDIR /Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/libs)
find_library(thermo_solver NAMES thermo_solver HINTS ${LIBSDIR})
find_library(superdrop_solver NAMES superdrop_solver HINTS ${LIBSDIR})
find_library(initialisation NAMES initialisation HINTS ${LIBSDIR})
target_link_libraries("${PROJECT_NAME}" PUBLIC thermo_solver superdrop_solver initialisation runmodel)

# # add directories to include for targets to us
target_include_directories(${PROJECT_NAME} PUBLIC ${LIBSDIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS OFF
  CXX_STANDARD 20)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -g)
#target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -pedantic -g)