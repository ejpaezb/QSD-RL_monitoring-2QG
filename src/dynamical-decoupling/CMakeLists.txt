cmake_minimum_required(VERSION 3.6)
project(IonTrap++)

# Notes:
# https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html

set(CMAKE_CXX_STANDARD 17) # enable C++20 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3 -flto")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -L /opt/homebrew/lib")

#SET(CMAKE_C_COMPILER mpicc)
#SET(CMAKE_CXX_COMPILER mpicxx)

# Bulding the external libraries
include(${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt)

include_directories(./experiment)
file(GLOB EXP_SOURCE_FILES "./experiment/*.h" "./experiment/*.cpp")

# ion-trap-qsd-1.3.5
#add_executable(IonTrap-QSD ion-trap-qsd-1.3.5/main.cpp ${EXP_SOURCE_FILES})
#target_link_libraries(IonTrap-QSD qsd cnpy boost_filesystem boost_serialization Eigen3::Eigen)

# ion-trap-odeint
#add_executable(IonTrap-ODE ion-trap-odeint/main.cpp ${EXP_SOURCE_FILES})
#target_link_libraries(IonTrap-ODE qsd cnpy boost_serialization)


# ion-trap-rl
find_package (Python3 COMPONENTS Interpreter Development NumPy)

#add_executable(IonTrap-COOLING ../cooling-dynamics/ion-trap-cooling/main.cpp ${EXP_SOURCE_FILES})
#target_link_libraries(IonTrap-COOLING qsd cnpy boost_filesystem boost_serialization Eigen3::Eigen)
#
add_executable(IonTrap-DD ion-trap-dyn_decoupling/main.cpp ${EXP_SOURCE_FILES})
target_link_libraries(IonTrap-DD qsd cnpy boost_filesystem boost_serialization Eigen3::Eigen)
#
#find_package(pybind11)
#include_directories(pybind11)
#
#add_executable(IonTrap-RL ion-trap-rl/main.cpp ${EXP_SOURCE_FILES})
#target_link_libraries(IonTrap-RL qsd cnpy boost_filesystem boost_serialization Eigen3::Eigen pybind11::embed)
