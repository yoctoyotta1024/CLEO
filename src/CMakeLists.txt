# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build
# or CXX=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# set project name and source directory of .cpp files
project("CLEOsetups")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# executable for cvodesdm coupled model
set(CVODECLEO "coupledCVODECLEO")
add_executable(${CVODECLEO} "main.cpp")
target_link_libraries(${CVODECLEO} PRIVATE
  cvodesdm superdrop_solver initialisation sdmgridboxes observers)
target_include_directories(${CVODECLEO} PRIVATE ${PROJECT_SOURCE_DIR})

# executable for cvodesdm coupled model
set(uncoupledCLEO "uncoupledCLEO")
add_executable(${uncoupledCLEO} "mainU.cpp")
target_link_libraries(${uncoupledCLEO} PRIVATE
  superdrop_solver initialisation sdmgridboxes observers)
target_include_directories(${uncoupledCLEO} PRIVATE ${PROJECT_SOURCE_DIR})

# set compiler properties for target(s)
set_target_properties(${CVODECLEO} ${uncoupledCLEO} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

target_compile_options(${CVODECLEO} PRIVATE -Wall -pedantic -g -gdwarf-4)
target_compile_options(${uncoupledCLEO} PRIVATE -Wall -pedantic -g -gdwarf-4)

#target_compile_options(${CVODECLEO} PRIVATE -Wall -Werror -pedantic -g -gdwarf-4)
#target_compile_options(${uncoupledCLEO} PRIVATE -Wall -Werror -pedantic -g -gdwarf-4)