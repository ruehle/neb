cmake_minimum_required(VERSION 2.8)

project(cpp_neb_tests)

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release Coverage." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)
ADD_DEFINITIONS(-std=c++0x)
ADD_DEFINITIONS(-Wall)



#cmake_policy(SET CMP0015 NEW)


# set the cpp_neb include directory
set(neb_include ../../src/neb)
include_directories(${neb_include})
message("neb include directory ${neb_include}")

# build the neb library
FILE(GLOB neb_sources ${neb_include}/*.cpp)
add_library(neb_lib ${neb_sources})

# get all the source files
FILE(GLOB sources *.cpp)


add_executable(run run_lj6.cpp)
target_link_libraries(run neb_lib)
