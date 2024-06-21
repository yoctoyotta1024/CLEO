# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set library name and print directory of this CMakeLists.txt (source directory of library) and print directory of this CMakeLists.txt (source directory of library)
set(LIBNAME "coupldyn_yac")
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIB_BINARY_DIR ${CLEO_BINARY_DIR}/lib)


message(STATUS "Using YAC installation from: ${YAC_ROOT}")
message(STATUS "Using YAXT installation from: ${YAXT_ROOT}")

find_package(MPI REQUIRED)
find_package(YAC REQUIRED)

# Add executables and create library target
set(SOURCES
"yac_comms.cpp"
"yac_cartesian_dynamics.cpp"
)
# must use STATIC (not(!) SHARED) lib for linking to executable if build is CUDA enabled with Kokkos
add_library("${LIBNAME}" STATIC ${SOURCES})
enable_language(C)

# Add directories for target library
# target_include_directories(${LIBNAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PRIVATE ${YAC_C_INCLUDE_DIR}) # YAC source directory
target_include_directories(${LIBNAME} PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# Link libraries to target library
set(LINKLIBS initialise gridboxes superdrops)
target_link_libraries("${LIBNAME}" PUBLIC ${LINKLIBS})
target_link_libraries("${LIBNAME}" PRIVATE YAC::YAC MPI::MPI_CXX)

# set specific C++ compiler options for target (optional)
#target_compile_options(${LIBNAME} PRIVATE)

# set C++ properties for target
set_target_properties(${LIBNAME} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
