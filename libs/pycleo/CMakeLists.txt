# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set library name and print directory of this CMakeLists.txt (source directory of library)
set(LIBNAME "pycleo")
message(STATUS "CLEO using ${LIBNAME} with LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIB_BINARY_DIR ${CLEO_BINARY_DIR}/pycleo)

# fetch and make pybind11 library available
include(FetchContent)
FetchContent_Declare(
  pybind11
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.13.6.tar.gz
  GIT_TAG a2e59f0e7065404b44dfe92a28aca47ba1378dc4
)
FetchContent_MakeAvailable(pybind11)

# # Include pybind11's CMake
# message(STATUS "pybind source dir: ${pybind11_SOURCE_DIR}")
# include(${pybind11_SOURCE_DIR}/CMake/pybind11.cmake)

message(STATUS "pybind11 installation in: ${CMAKE_CURRENT_BINARY_DIR}")

# Add executables and create library target
set(SOURCES
"pycleo.cpp"
)
pybind11_add_module(${LIBNAME} MODULE ${SOURCES})

# Add directories for target library
#target_include_directories(${LIBNAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# Link libraries to target library
set(LINKLIBS gridboxes superdrops runcleo)
target_link_libraries(${LIBNAME} PUBLIC "${LINKLIBS}")
target_link_libraries(${LIBNAME} PUBLIC Kokkos::kokkos)

# set specific C++ compiler options for target (optional)
#target_compile_options(${LIBNAME} PRIVATE)

# set C++ properties for target
set_target_properties(${LIBNAME} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  ARCHIVE_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
