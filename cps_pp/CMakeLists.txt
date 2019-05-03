cmake_minimum_required(VERSION 2.8)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
set(CMAKE_COLOR_MAKEFILE ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(include)

find_package(MPI)

# include_directories(SYSTEM "/home/jiquntu/cps-build/cmake/local/include")
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

set( ext_loc ${CMAKE_BINARY_DIR}/deps )

include(ExternalProject)

ExternalProject_Add( gsl
  URL http://gnu.mirrors.pair.com/gsl/gsl-1.16.tar.gz
  PREFIX ${CMAKE_BINARY_DIR}/deps/gsl
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${ext_loc}
  BUILD_COMMAND make -j60
  INSTALL_COMMAND make install
)

ExternalProject_Add( qmp
  URL https://github.com/usqcd-software/qmp/archive/qmp2-3-1-CJ.zip
  PREFIX ${CMAKE_BINARY_DIR}/deps/qmp
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${ext_loc} --with-qmp-comms-type=mpi CXX=mpicxx CC=mpicc
  BUILD_COMMAND make
  INSTALL_COMMAND make install
)

ExternalProject_Add( qio
  URL http://usqcd-software.github.io/downloads/qio/qio-2.3.8.tar.gz
  PREFIX ${CMAKE_BINARY_DIR}/deps/qio
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${ext_loc} --build=none --with-qmp=${ext_loc} --enable-largefile CFLAGS=-O2
  BUILD_COMMAND make -j60
  INSTALL_COMMAND make install
)

add_dependencies( qio qmp )

include_directories(${ext_loc}/include)

# The arg files
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/generated/include/alg)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/generated/src)
include_directories(${CMAKE_BINARY_DIR}/generated/include)
add_subdirectory(misc)

add_definitions(-g -O2 -fopenmp)
add_definitions(-DUSE_OMP -DUSE_QMP -DUSE_QIO)

# The main library
add_subdirectory(src)

