#
# CMake for EigenLab
#

cmake_minimum_required (VERSION 2.8)

project (EigenLab CXX)

set (EIGEN_VERSION "3.3.4" CACHE STRING "Select Eigen Version.")

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

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-int-in-bool-context") # eigen 3.3.4
add_definitions (-DDEBUG)
add_definitions (-D_GLIBCXX_DEBUG)
add_executable (Test Test.cpp)
add_dependencies (Test eigenX)

#
# Afl
#
ExternalProject_Add (aflX
  PREFIX aflX
  URL http://lcamtuf.coredump.cx/afl/releases/afl-2.52b.tgz
  CONFIGURE_COMMAND ""
  BUILD_COMMAND "make"
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  )
ExternalProject_Get_Property (aflX source_dir)
set (AFL_DIR "${source_dir}")
message (${AFL_DIR})

#
# Any add_executable after this point will be compiled with AFL tooling!
#
set (CMAKE_CXX_COMPILER ${AFL_DIR}/afl-g++)

add_executable (EigenLab EigenLab.cpp)
add_dependencies (EigenLab eigenX)

add_custom_target (fuzz
  COMMAND ${AFL_DIR}/afl-fuzz -i fuzz-data -o fuzz-out $<TARGET_FILE:EigenLab> -f @@
  )


enable_testing ()
add_test (basic Test)

#add_executable (tmp tmp.cpp)
