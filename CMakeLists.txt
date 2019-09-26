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
find_package (Threads)
link_libraries (${CMAKE_THREAD_LIBS_INIT})

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

# TBB
#find_package(TBB REQUIRED)
#add_definitions(${TBB_DEFINITIONS})
#include_directories(${TBB_INCLUDE_DIRS})
#link_libraries(${TBB_LIBRARIES})

# MKL
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIR})
link_libraries(${MKL_LIBRARIES})

# BLAS Math Libraries
find_package(BLAS REQUIRED)
link_libraries(${BLAS_LIBRARIES})
  
find_package(LAPACK REQUIRED)
link_libraries(${LAPACK_LIBRARIES})

# Include directories
include_directories(
  ${pipe_network_SOURCE_DIR}/include/
  ${pipe_network_SOURCE_DIR}/external/
)

link_directories(${pipe_network_SOURCE_DIR}/external/)

# pipe_network executable
SET(pipe_network_src
  src/curves.cc
  src/hydralic_sim.cc
  src/input.cc
  src/matrix_assembler.cc
  src/mesh.cc
  src/mkl_unsym.cc
  src/pardiso_unsym.cc
  src/output.cc
  )


add_library(lpipe_network SHARED ${pipe_network_src})

add_executable(pipe_network ${pipe_network_SOURCE_DIR}/src/main.cc)
target_link_libraries(pipe_network lpipe_network)

if(PIPE_NETWORK_TEST_BUILD_TESTING)
  # Unit test
  SET(pipe_network_test_src
    tests/mesh_test.cc
    tests/node_test.cc
    tests/link_test.cc
    tests/matrix_assembler_test.cc
    tests/output_test.cc
    tests/hydraulic_sim_test.cc
    tests/input_test.cc
    tests/curve_test.cc
    )
  add_executable(pipe_network_test tests/main_test.cc ${pipe_network_test_src})
  target_link_libraries(pipe_network_test lpipe_network)
  add_test(NAME pipe_network_test COMMAND $<TARGET_FILE:pipe_network_test>)
  enable_testing()
endif()

# Coverage
find_package(codecov)
if(ENABLE_COVERAGE)
  add_executable(pipe_network_test_coverage tests/main_test.cc ${pipe_network_src} ${pipe_network_test_src})
  add_coverage(pipe_network_test_coverage)
endif()
