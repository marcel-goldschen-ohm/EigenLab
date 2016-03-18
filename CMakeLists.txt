#
# CMake for EigenLab
#

cmake_minimum_required (VERSION 2.8)
set (Boost_NO_BOOST_CMAKE ON) # hack for stupid old broken cmakes

project (EigenLab CXX)

find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
#find_package (Eigen3 REQUIRED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -std=c++11")

include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${EIGEN3_INCLUDE_DIRS}")

add_executable (Test Test.cpp)
target_compile_definitions (Test PRIVATE DEBUG)

enable_testing ()
add_test (basic Test)

#add_executable (tmp tmp.cpp)
