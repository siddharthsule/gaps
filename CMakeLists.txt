# Minimum required version of CMake
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# Project name and languages used
project(GPUEvGen LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define directories
set(BINDIR ${CMAKE_SOURCE_DIR}/bin)
set(INCDIR ${CMAKE_SOURCE_DIR}/include)
set(SRCDIR ${CMAKE_SOURCE_DIR}/src)

# Define source and include files
file(GLOB SRCS "${SRCDIR}/*.cpp" "${SRCDIR}/*.cu")
file(GLOB INCS "${INCDIR}/*.h" "${INCDIR}/*.cuh")

# Set all source files to compile with CUDA
set_source_files_properties(${SRCS} PROPERTIES LANGUAGE CUDA)

# Include directories
include_directories(${INCDIR})

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINDIR})

# Add executable
add_executable(GPUEvGen ${SRCS})

# Set CUDA architecture off, so it uses the one available
set_property(TARGET GPUEvGen PROPERTY CUDA_ARCHITECTURES OFF)

# Set CUDA properties
set_target_properties(
    GPUEvGen
    PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)