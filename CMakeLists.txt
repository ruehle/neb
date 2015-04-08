cmake_minimum_required(VERSION 2.6)

project(neb)

set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "ruehle@vota.org")

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)
ADD_DEFINITIONS(-Wall)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)


