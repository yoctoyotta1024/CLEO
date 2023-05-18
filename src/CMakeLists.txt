# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-13 CC=/opt/homebrew/bin/gcc-13 cmake -S ./ -B ./build
# or CXX=g++ CC=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) # for using Kokkos with NVC++ compiler 

# set project name and source directory of .cpp files
project("CLEOsetups")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Kokkos library stuff
set(CMAKE_CXX_FLAGS_B4 "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD_B4 "${CMAKE_CXX_STANDARD}")
include(FetchContent)
FetchContent_Declare(
  Kokkos
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/kokkos/kokkos/archive/refs/tags/4.0.01.tar.gz
  GIT_TAG 62d2b6c879b74b6ae7bd06eb3e5e80139c4708e6
)
set(CMAKE_CXX_STANDARD "20")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
FetchContent_MakeAvailable(Kokkos)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_B4}")
set(CMAKE_CXX_STANDARD "${CMAKE_CXX_STANDARD_B4}")

set(SDMLIBS superdrop_solver initialisation sdmgridboxes observers)
#set(compileops -Werror -Wall -pedantic -g -gdwarf-4)
set(compileops -Wall -pedantic -g -gdwarf-4) # less strict without -Werror
#set(compileops  -Wall -pedantic -g) # cannot use -gdwarf-4 flag for nvc++ compiler 

# primary executable for CLEO sdm
set(runCLEO "runCLEO")
add_executable(${runCLEO} "main.cpp")
target_link_libraries(${runCLEO} PRIVATE thermofromfile "${SDMLIBS}")
target_link_libraries(${runCLEO} PUBLIC Kokkos::kokkos)
target_include_directories(${runCLEO} PRIVATE ${PROJECT_SOURCE_DIR})
target_compile_options(${runCLEO} PRIVATE ${compileops})

# executable for CVODE-CLEO coupled model 0D golovin collisions setup
set(golcolls0D "golcolls0D")
add_executable(${golcolls0D} EXCLUDE_FROM_ALL "main_golcolls0D.cpp")
target_link_libraries(${golcolls0D} PRIVATE thermofromfile "${SDMLIBS}")
target_link_libraries(${golcolls0D} PUBLIC Kokkos::kokkos)
target_compile_options(${golcolls0D} PRIVATE ${compileops})

# executable for CVODE-CLEO coupled model 0D condensation setup 
set(cond0D "cond0D")
add_executable(${cond0D} EXCLUDE_FROM_ALL "main_cond0D.cpp")
target_link_libraries(${cond0D} PRIVATE cvodecoupld "${SDMLIBS}")
target_link_libraries(${cond0D} PUBLIC Kokkos::kokkos)
target_compile_options(${cond0D} PRIVATE ${compileops}) 

# set compiler properties for target(s)
set_target_properties(${runCLEO} ${cond0D} ${golcolls0D} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)