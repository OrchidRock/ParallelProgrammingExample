add_executable(${PROJECT_NAME}_serial serial_main.cpp)
target_link_libraries(${PROJECT_NAME}_serial PRIVATE gmp gmpxx)

add_executable(${PROJECT_NAME}_openmp openmp_main.cpp)
add_compile_options(-fopenmp)
target_link_libraries(${PROJECT_NAME}_openmp PRIVATE gmp gmpxx)
