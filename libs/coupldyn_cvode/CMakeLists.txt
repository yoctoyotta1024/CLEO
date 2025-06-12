# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set library name and print directory of this CMakeLists.txt (source directory of library) and print directory of this CMakeLists.txt (source directory of library)
set(LIBNAME "coupldyn_cvode")
message(STATUS "CLEO using ${LIBNAME} with LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIB_BINARY_DIR ${CLEO_BINARY_DIR}/lib)

# fetch CVODE libraries and make available
include(FetchContent)
FetchContent_Declare(
  cvodes
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/LLNL/sundials/releases/download/v6.5.0/cvodes-6.5.0.tar.gz
  URL_HASH MD5=dc0cd29340599f0f8b027c6f1e0107ba
)
FetchContent_MakeAvailable(cvodes)

# Add executables and create library target
set(SOURCES
"cvodecomms.cpp"
"cvodedynamics.cpp"
"differentialfuncs.cpp"
)
# must use STATIC (not(!) SHARED) lib for linking to executable if build is CUDA enabled with Kokkos
add_library("${LIBNAME}" STATIC ${SOURCES})

# Add directories for target library
#target_include_directories(${LIBNAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# Link libraries to target library
set(LINKLIBS configuration superdrops gridboxes)
target_link_libraries("${LIBNAME}" PUBLIC ${LINKLIBS})
target_link_libraries("${LIBNAME}" PUBLIC SUNDIALS::cvodes_static Kokkos::kokkos)

# set specific C++ compiler options for target (optional)
#target_compile_options(${LIBNAME} PRIVATE)

# set C++ properties for target
set_target_properties(${LIBNAME} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  ARCHIVE_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
