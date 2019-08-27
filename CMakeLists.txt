project(pipe_network LANGUAGES CXX)

# Require C++14-compliant compiler; only available for CMake v. 3.1 and up
set(CMAKE_CXX_STANDARD 14)

cmake_minimum_required(VERSION 3.1)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 SET(CMAKE_BUILD_TYPE "Debug")
 #SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# GNU Specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -pthread")
endif()

# Clang Specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -pthread")
endif()

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Boost Archive
#find_package(Boost REQUIRED COMPONENTS system serialization)
#include_directories(${BOOST_INCLUDE_DIRS})
#link_libraries(${Boost_SERIALIZATION_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Pthreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# TBB
#find_package(TBB REQUIRED)
#add_definitions(${TBB_DEFINITIONS})
#include_directories(${TBB_INCLUDE_DIRS})
#link_libraries(${TBB_LIBRARIES})

# Include directories
include_directories(
  ${pipe_network_SOURCE_DIR}/include/
  ${pipe_network_SOURCE_DIR}/external/
)

link_directories(${pipe_network_SOURCE_DIR}/external/)

# pipe_network executable
SET(pipe_network_src
#        src/main.cc
#        src/mesh.cc
#        src/matrix_assembler.cc
#        src/eigen_cg.cc
#        src/hydralic_sim.cc
#        src/eigen_gmres.cc
#        src/pardiso_unsym.cc
        src/pardiso_unsym.cpp)
#        src/input.cc)
add_executable(pipe_network
  ${pipe_network_src}
)

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(pipe_network ${lapackblas_libraries})
target_link_libraries(pipe_network OpenMP::OpenMP_CXX)
target_link_libraries(pipe_network Threads::Threads)
target_link_libraries(pipe_network libpardiso600-MACOS-X86-64.dylib)

# Unit test
SET(test_src
        src/mesh.cc
        src/matrix_assembler.cc
        src/hydralic_sim.cc
        src/eigen_gmres.cc
        src/input.cc
        src/eigen_cg.cc
        src/eigen_cg_lem.cc
        src/pardiso_unsym.cc
        src/curves.cc
        src/output.cc
        tests/mesh_test.cc
        tests/node_test.cc
        tests/link_test.cc
        tests/matrix_assembler_test.cc
#        tests/output_test.cc
        tests/hydraulic_sim_test.cc
        tests/input_test.cc
#        tests/curve_test.cc
        tests/solver_test.cc
        tests/test.cc
        )
add_executable(pipe_network_test
  ${test_src}
)

target_link_libraries(pipe_network_test ${lapackblas_libraries})
target_link_libraries(pipe_network_test OpenMP::OpenMP_CXX)
target_link_libraries(pipe_network_test Threads::Threads)
target_link_libraries(pipe_network_test libpardiso600-MACOS-X86-64.dylib)
target_link_libraries(pipe_network_test libgfortran.dylib)

add_test(NAME pipe_network_test COMMAND $<TARGET_FILE:pipe_network_test>)
enable_testing()

# Coverage
find_package(codecov)
add_coverage(pipe_network)
add_coverage(pipe_network_test)
