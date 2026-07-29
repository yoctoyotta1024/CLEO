# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# fetch and make yaml-cpp library available
include(FetchContent)
FetchContent_Declare(
  yaml-cpp
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.9.0.tar.gz
  GIT_TAG 56e3bb550c91fd7005566f19c079cb7a503223cf
)
FetchContent_MakeAvailable(yaml-cpp)

message(STATUS "yaml-cpp installation in: ${CMAKE_CURRENT_BINARY_DIR}")
