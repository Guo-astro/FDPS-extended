cmake_minimum_required(VERSION 3.15)
project(FDPS_Tests CXX)

# Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# FDPS source directory
set(FDPS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)

# Include FDPS headers
include_directories(${FDPS_SRC_DIR})

# Google Test setup
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Enable testing
enable_testing()
include(GoogleTest)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Unit System and Output Tests
add_executable(unit_system_gtest
    unit_system_gtest.cpp
)

target_link_libraries(unit_system_gtest
    GTest::gtest_main
    GTest::gmock
)

# Optional: HDF5 support
find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
    target_compile_definitions(unit_system_gtest PRIVATE PARTICLE_SIMULATOR_USE_HDF5)
    target_include_directories(unit_system_gtest PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(unit_system_gtest ${HDF5_LIBRARIES})
    message(STATUS "HDF5 found - enabling HDF5 tests")
else()
    message(STATUS "HDF5 not found - HDF5 tests will be skipped")
endif()

# Discover tests
gtest_discover_tests(unit_system_gtest
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    PROPERTIES LABELS "unit"
)

# Add custom target to run all tests
add_custom_target(check
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --verbose
    DEPENDS unit_system_gtest
    COMMENT "Running all tests"
)

# Install test executable (optional)
install(TARGETS unit_system_gtest
    RUNTIME DESTINATION bin/tests
)
