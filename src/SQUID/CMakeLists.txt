cmake_minimum_required(VERSION 3.6)
project(CVCQED++)

# Notes:
# 1. Set the proper -mtune and -march arguments when running on computecanada
# 2. Set the -g argument for debugging purposes

set(CMAKE_CXX_STANDARD 11) # enable C++11 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast") # https://wiki.gentoo.org/wiki/GCC_optimization
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L /usr/local/lib")

# Bulding the external libraries
include(${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt)

include_directories(./experiment/)
file(GLOB EXP_SOURCE_FILES "./experiment/*.h" "./experiment/*.cpp")

# SQUID-qsd-1.3.5
add_executable(CV-QSD cv-qsd-1.3.5/main.cpp ${EXP_SOURCE_FILES})
target_link_libraries(CV-QSD qsd cnpy boost_serialization)
