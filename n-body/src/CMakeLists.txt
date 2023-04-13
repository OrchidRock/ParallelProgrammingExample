add_executable(${PROJECT_NAME}_serial serial_main.cpp)
target_link_libraries(${PROJECT_NAME}_serial PRIVATE gmp gmpxx)

add_executable(${PROJECT_NAME}_openmp openmp_main.cpp)
target_include_directories(${PROJECT_NAME}_openmp PRIVATE ${OpenMP_CXX_INCLUDE_DIRS})
target_compile_options(${PROJECT_NAME}_openmp PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME}_openmp PRIVATE gmp gmpxx ${OpenMP_CXX_LIBRARIES})


add_executable(${PROJECT_NAME}_mpi mpi_main.cpp)
target_include_directories(${PROJECT_NAME}_mpi PRIVATE ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_mpi PRIVATE gmp gmpxx ${MPI_CXX_LIBRARIES})

add_executable(${PROJECT_NAME}_cuda cuda_main.cu)
target_link_libraries(${PROJECT_NAME}_cuda PRIVATE gmp gmpxx)
