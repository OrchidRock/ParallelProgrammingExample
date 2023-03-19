add_executable(${PROJECT_NAME}_serial serial_main.cpp)
target_link_libraries(${PROJECT_NAME}_serial PRIVATE gmp gmpxx)

add_executable(${PROJECT_NAME}_openmp openmp_main.cpp)
add_compile_options(-fopenmp)
target_link_libraries(${PROJECT_NAME}_openmp PRIVATE gmp gmpxx)

add_executable(${PROJECT_NAME}_mpi mpi_main.cpp)
target_include_directories(${PROJECT_NAME}_mpi PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_mpi PRIVATE gmp gmpxx ${MPI_CXX_LIBRARIES})
