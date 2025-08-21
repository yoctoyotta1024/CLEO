# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  # cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("boxmodelcollisions")
message(STATUS "CLEO including ${PROJECT_NAME} with PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)

# create primary executable for CLEO uncoupled model golovin kernel box model setup
add_executable(golcolls EXCLUDE_FROM_ALL "main_golcolls.cpp")
add_executable(longcolls EXCLUDE_FROM_ALL "main_longcolls.cpp")
add_executable(lowlistcolls EXCLUDE_FROM_ALL "main_lowlistcolls.cpp")
add_executable(szakallurbichcolls EXCLUDE_FROM_ALL "main_szakallurbichcolls.cpp")
add_executable(testikstraubcolls EXCLUDE_FROM_ALL "main_testikstraubcolls.cpp")

# Add directories and link libraries to target
target_link_libraries(golcolls PRIVATE cartesiandomain "${CLEOLIBS}")
target_link_libraries(longcolls PRIVATE cartesiandomain "${CLEOLIBS}")
target_link_libraries(lowlistcolls PRIVATE cartesiandomain "${CLEOLIBS}")
target_link_libraries(szakallurbichcolls PRIVATE cartesiandomain "${CLEOLIBS}")
target_link_libraries(testikstraubcolls PRIVATE cartesiandomain "${CLEOLIBS}")

target_link_libraries(golcolls PUBLIC Kokkos::kokkos)
target_link_libraries(longcolls PUBLIC Kokkos::kokkos)
target_link_libraries(lowlistcolls PUBLIC Kokkos::kokkos)
target_link_libraries(szakallurbichcolls PUBLIC Kokkos::kokkos)
target_link_libraries(testikstraubcolls PUBLIC Kokkos::kokkos)

target_include_directories(golcolls PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory
target_include_directories(longcolls PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory
target_include_directories(lowlistcolls PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory
target_include_directories(szakallurbichcolls PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory
target_include_directories(testikstraubcolls PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# set specific C++ compiler options for target (optional)
#target_compile_options([target] PRIVATE)

# set compiler properties for target(s)
set_target_properties(golcolls PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
set_target_properties(longcolls PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
set_target_properties(lowlistcolls PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
set_target_properties(szakallurbichcolls PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
set_target_properties(testikstraubcolls PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
