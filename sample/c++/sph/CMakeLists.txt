cmake_minimum_required(VERSION 3.15)
project(FDPS_SPH CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add FDPS include directory
include_directories(../../../src/)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -funroll-loops")

# Check for OpenMP support
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_compile_definitions(PARTICLE_SIMULATOR_THREAD_PARALLEL)
    message(STATUS "OpenMP found and enabled")
else()
    message(STATUS "OpenMP not found, building without thread parallelization")
endif()

# Check for MPI support
find_package(MPI)
if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    add_compile_definitions(PARTICLE_SIMULATOR_MPI_PARALLEL)
    message(STATUS "MPI found and enabled")
else()
    message(STATUS "MPI not found, building without MPI parallelization")
endif()

# Source files
set(SOURCES main.cpp)

# Create executable
add_executable(sph.out ${SOURCES})

# Link libraries
if(MPI_CXX_FOUND)
    target_link_libraries(sph.out ${MPI_CXX_LIBRARIES})
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(sph.out OpenMP::OpenMP_CXX)
endif()
