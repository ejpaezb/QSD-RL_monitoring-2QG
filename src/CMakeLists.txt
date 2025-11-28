# Bulding the external libraries
find_package(Boost)
IF (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
endif()

# set(Eigen3_DIR path_to_eigen)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(../../include/)

file(GLOB QSD_SOURCE_FILES "../../include/qsd/*.h" "../../include/qsd/*.cpp")
file(GLOB CNPY_SOURCE_FILES "../../include/cnpy/*.h" "../../include/cnpy/*.cpp")


add_library(qsd STATIC ${QSD_SOURCE_FILES})
add_library(cnpy STATIC ${CNPY_SOURCE_FILES})

#add_library(xyz INTERFACE)
#target_include_directories(de INTERFACE path-to-xyz)