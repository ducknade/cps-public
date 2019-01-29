cmake_minimum_required(VERSION 2.8)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
set(CMAKE_COLOR_MAKEFILE ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(include)

find_package(MPI)

add_definitions(-g -O2 -fopenmp)
add_definitions(-DUSE_OMP -DUSE_QMP -DUSE_QIO)

include_directories(SYSTEM "/home/jiquntu/cps-build/cmake/local/include")
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

add_subdirectory(misc)

add_subdirectory(src)
