project(pipe_network LANGUAGES CXX)

# Require C++17-compliant compiler; only available for CMake v. 3.1 and up
set(CMAKE_CXX_STANDARD 17)

cmake_minimum_required(VERSION 3.1)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 SET(CMAKE_BUILD_TYPE "Debug")
# SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# GNU Specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -pthread")
endif()

# Clang Specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -pthread")
endif()

# CMake seems to have no way to enable/disable testing per subproject,
# so we provide an option similar to BUILD_TESTING, but just for CITYSCAPE.
option(PIPE_NETWORK_TEST_BUILD_TESTING "enable testing for cityscape" ON)

SET(USE_PARDISO  OFF)
SET(USE_CUDA  OFF)

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# pipe_network executable
SET(pipe_network_src
        src/index_manager.cc
        src/mesh_components.cc
        src/mesh.cc
        src/matrix_assembler_components.cc
        src/matrix_assembler.cc
        src/curves.cc

        )
#        src/curves.cc
#        src/hydralic_sim.cc
#        src/input.cc
#        src/matrix_assembler.cc
#        src/mesh.cc
#        src/output.cc
#        src/mkl_unsym.cc
#        src/pardiso_unsym.cc)




# Find dependencies for the library
# Eigen
find_package(Eigen3 REQUIRED)

# Boost Archive
find_package(Boost COMPONENTS system filesystem REQUIRED)


# Pthreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads)

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

# MKL
find_package(MKL REQUIRED)

# BLAS Math Libraries
find_package(BLAS REQUIRED)

# LAPACK
find_package(LAPACK REQUIRED)


# Solver options
if(USE_CUDA)
  find_package(CUDA REQUIRED)
  if (CUDA_FOUND)
    SET(CUDA_LIB_USED ${CUDA_CUDART_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cusparse_LIBRARY})
    SET(NVIDIA_HELPERS "/home/rewu/Documents/research/NVIDIA_CUDA-10.0_Samples/common/inc")
    # create library
    add_library(lpipe_network SHARED ${pipe_network_src} src/cuda_unsym.cc)
    target_include_directories(lpipe_network
            PUBLIC
            ${NVIDIA_HELPERS}
            ${CUDA_INCLUDE_DIRS})
    target_link_libraries(lpipe_network PUBLIC  ${CUDA_LIB_USED})
  endif ()
  else()
  add_library(lpipe_network SHARED ${pipe_network_src})
endif()

# PARDISO6
if(USE_PARDISO)
  SET(PARDISO6 "${pipe_network_SOURCE_DIR}/external/libpardiso600-GNU720-X86-64.so")
  target_link_libraries(lpipe_network PUBLIC ${PARDISO6})
  else()
  target_link_libraries(lpipe_network PUBLIC ${MKL_LIBRARIES})
endif()

# include directories
target_include_directories(lpipe_network
        PUBLIC
        ${pipe_network_SOURCE_DIR}/include/
        ${pipe_network_SOURCE_DIR}/external/
        ${EIGEN3_INCLUDE_DIR}
        ${MKL_INCLUDE_DIR}
        ${BOOST_INCLUDE_DIRS}

        )

# link libs
target_link_libraries(lpipe_network
        PUBLIC
        ${CMAKE_THREAD_LIBS_INIT}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${Boost_LIBRARIES}
        )

#add_executable(pipe_network ${pipe_network_SOURCE_DIR}/src/main.cc)
#target_link_libraries(pipe_network
#        PUBLIC
#        lpipe_network)

# Testing
if(PIPE_NETWORK_TEST_BUILD_TESTING)
  # Unit test
  SET(pipe_network_test_src
          tests/node_test.cc
          tests/link_test.cc
          tests/mesh_test.cc
          tests/matrix_assembler_test.cc
#          tests/solver_test.cc
#          tests/input_test.cc
#          tests/curve_test.cc
#          tests/hydraulic_sim_test.cc
#          tests/sim_pressure_test.cc
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
  target_link_libraries(pipe_network_test_coverage lpipe_network)
  add_coverage(pipe_network_test_coverage)
endif()
