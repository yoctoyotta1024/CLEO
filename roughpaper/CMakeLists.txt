# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# e.g. could run cmake with these compiler options:
# CXX=g++-13 CC=gcc-13 cmake -S ./ -B ./build -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON
# CXX=g++ CC=gcc cmake -S ./ -B ./build -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON

# set cmake version
cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) # for using Kokkos with NVC++ compiler 

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
include_directories("${CLEODIR}/libs")
# include_directories(${PROJECT_SOURCE_DIR})
# include_directories(${CLEODIR})

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/initialisation" ${CMAKE_CURRENT_BINARY_DIR}/libs/initialisation)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/superdrop_solver" ${CMAKE_CURRENT_BINARY_DIR}/libs/superdrop_solver)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/observers" ${CMAKE_CURRENT_BINARY_DIR}/libs/observers)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/sdmgridboxes" ${CMAKE_CURRENT_BINARY_DIR}/libs/sdmgridboxes)

set(LINKLIBS superdrop_solver initialisation sdmgridboxes observers)
target_link_libraries("${PROJECT_NAME}" PUBLIC ${LINKLIBS})

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

# link kokkos library to target
target_link_libraries(${PROJECT_NAME} PUBLIC Kokkos::kokkos)

set_target_properties(${PROJECT_NAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

#set(compileops -Werror -Wall -pedantic -g -gdwarf-4)
set(compileops -Wall -pedantic -g -gdwarf-4) # less strict without -Werror
#set(compileops  -Wall -pedantic -g) # cannot use -gdwarf-4 flag for nvc++ compiler 

target_compile_options(${PROJECT_NAME} PUBLIC ${compileops})