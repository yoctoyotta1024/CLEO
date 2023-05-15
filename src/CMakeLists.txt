# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-13 CC=/opt/homebrew/bin/gcc-13 cmake -S ./ -B ./build
# or CXX=g++ CC=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# set project name and source directory of .cpp files
project("CLEOsetups")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

set(SDMLIBS superdrop_solver initialisation sdmgridboxes observers)
set(compileops -Wall -pedantic -g -gdwarf-4)
#set(compileops -Werror -Wall -pedantic -g -gdwarf-4)

# primary executable for CLEO sdm
set(runCLEO "runCLEO")
add_executable(${runCLEO} "main.cpp")
target_link_libraries(${runCLEO} PRIVATE thermofromfile "${SDMLIBS}")
target_include_directories(${runCLEO} PRIVATE ${PROJECT_SOURCE_DIR})
target_compile_options(${runCLEO} PRIVATE ${compileops})

# executable for CVODE-CLEO coupled model 0D golovin collisions setup
set(golcolls0D "golcolls0D")
add_executable(${golcolls0D} EXCLUDE_FROM_ALL "main_golcolls0D.cpp")
target_link_libraries(${golcolls0D} PRIVATE thermofromfile "${SDMLIBS}")
target_compile_options(${golcolls0D} PRIVATE ${compileops})

# executable for CVODE-CLEO coupled model 0D condensation setup 
set(cond0D "cond0D")
add_executable(${cond0D} EXCLUDE_FROM_ALL "main_cond0D.cpp")
target_link_libraries(${cond0D} PRIVATE cvodecoupld "${SDMLIBS}")
target_compile_options(${cond0D} PRIVATE ${compileops}) 


# set compiler properties for target(s)
set_target_properties(${runCLEO} ${cond0D} ${golcolls0D} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)