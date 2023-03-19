cmake_minimum_required(VERSION 3.5.1)
project(TSP VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(MPI REQUIRED)
#message(STATUS "${MPI_CXX_INCLUDE_DIRS} ${MPI_CXX_LIBRARIES} ${MPI_CXX_COMPILE_OPTIONS}")
add_subdirectory(src)
