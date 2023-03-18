add_executable(${PROJECT_NAME}_serial serial_main.cpp)
target_link_libraries(${PROJECT_NAME}_serial PRIVATE gmp gmpxx)
