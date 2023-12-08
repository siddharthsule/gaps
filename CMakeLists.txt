cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(GPUEvGen LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Directories
set(BINDIR ${CMAKE_SOURCE_DIR}/bin)
set(INCDIR ${CMAKE_SOURCE_DIR}/include)
set(SRCDIR ${CMAKE_SOURCE_DIR}/src)

# Files
file(GLOB SRCS "${SRCDIR}/*.cpp" "${SRCDIR}/*.cu")
file(GLOB INCS "${INCDIR}/*.h" "${INCDIR}/*.cuh")

# Set all source files to CUDA language
set_source_files_properties(${SRCS} PROPERTIES LANGUAGE CUDA)

# Include directories
include_directories(${INCDIR})

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINDIR})

# Add executable
add_executable(GPUEvGen ${SRCS})

# Link CUDA libraries
set_target_properties(
    GPUEvGen
    PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)