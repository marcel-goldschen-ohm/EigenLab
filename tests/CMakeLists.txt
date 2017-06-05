#
# CMake for EigenLab
#

cmake_minimum_required (VERSION 2.8)

project (EigenLab CXX)

set (EIGEN_VERSION "3.3.3" CACHE STRING "Select Eigen Version.")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

message ("Testing EigenLab against Eigen-${EIGEN_VERSION}")

#
# Eigen
#
include (ExternalProject)
ExternalProject_Add (eigenX
  PREFIX eigenX
  URL http://bitbucket.org/eigen/eigen/get/${EIGEN_VERSION}.tar.bz2
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property (eigenX source_dir)
set (EIGEN3_INCLUDE_DIRS "${source_dir}")

include_directories ("${PROJECT_SOURCE_DIR}/..")
include_directories ("${EIGEN3_INCLUDE_DIRS}")

add_definitions (-DDEBUG)
add_definitions (-D_GLIBCXX_DEBUG)
add_executable (Test Test.cpp)
add_dependencies (Test eigenX)

enable_testing ()
add_test (basic Test)

#add_executable (tmp tmp.cpp)
