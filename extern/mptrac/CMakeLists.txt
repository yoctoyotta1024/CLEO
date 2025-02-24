# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# fetch and make mptrac library available
include(FetchContent)
FetchContent_Declare(
  mptrac
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/slcs-jsc/mptrac/archive/refs/heads/ww_dd.zip
  # GIT_TAG 0bff3c8e82c4b6b9b48e6472898682d1d50e7f2b
)
FetchContent_MakeAvailable(mptrac)

message(STATUS "mptrac installation in: ${CMAKE_CURRENT_BINARY_DIR}")
