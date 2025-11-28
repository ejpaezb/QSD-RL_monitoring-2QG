cmake_minimum_required(VERSION 3.6)
project(QSD-Playground)

# Notes:
# 1. Set the proper -mtune and -march arguments when running on computecanada
# 2. Set the -g argument for debugging purposes

set(CMAKE_CXX_STANDARD 11) # enable C++11 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast") # https://wiki.gentoo.org/wiki/GCC_optimization
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L /usr/local/lib")

#SET(CMAKE_C_COMPILER mpicc)
#SET(CMAKE_CXX_COMPILER mpicxx)

# Bulding the external libraries
include(${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt)

# qsd-playground
add_executable(QSD-Playground main.cpp)
target_link_libraries(QSD-Playground qsd cnpy boost_serialization Eigen3::Eigen)