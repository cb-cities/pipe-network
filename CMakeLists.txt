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

# Include directories
include_directories(
  ${pipe_network_SOURCE_DIR}/include/
  ${pipe_network_SOURCE_DIR}/external/
)

# pipe_network executable
SET(pipe_network_src
  src/main.cc
  src/eigen_gmres.cc
  src/matrix_assembler.cc
  src/mesh.cc
  src/pipe.cc
        src/node.cc
  src/io.cc
)
add_executable(pipe_network
  ${pipe_network_src}
)

target_link_libraries(pipe_network Threads::Threads)

# Unit test
SET(test_src
   src/eigen_gmres.cc
   src/matrix_assembler.cc
   src/mesh.cc
   src/pipe.cc
   src/io.cc
        src/node.cc
#   tests/eigen_gmres_test.cc
#   tests/matrix_assembler_test.cc
#        tests/matrix_assembler_test_full.cc
#   tests/mesh_test.cc
#   tests/node_test.cc
#   tests/pipe_test.cc
#   tests/io_test.cc
   tests/todini_example_test.cc
#        tests/todini_no_demand.cc
   tests/test.cc
)   
add_executable(pipe_network_test
  ${test_src}
)

target_link_libraries(pipe_network_test Threads::Threads)

add_test(NAME pipe_network_test COMMAND $<TARGET_FILE:pipe_network_test>)
enable_testing()

# Coverage
find_package(codecov)
add_coverage(pipe_network)
add_coverage(pipe_network_test)
