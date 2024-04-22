cmake_minimum_required(VERSION 3.10)
project(base)

# Enable CUDA
enable_language(CUDA) 

set(CMAKE_CXX_STANDARD 17)

include_directories(include)
file(GLOB SOURCES "src/*.cu")

# Set all source files to compile with CUDA
set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CUDA)

add_library(base ${SOURCES})

# Set CUDA architecture to 7.0 for Tesla V100
set_property(TARGET base PROPERTY CUDA_ARCHITECTURES 70)