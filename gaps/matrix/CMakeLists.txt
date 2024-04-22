cmake_minimum_required(VERSION 3.10)
project(matrix)

# Enable CUDA
enable_language(CUDA) 

set(CMAKE_CXX_STANDARD 17)

include_directories(include ../base/include)
file(GLOB SOURCES "src/*.cu")

# Set all source files to compile with CUDA
set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CUDA)

add_library(matrix ${SOURCES})

# Set CUDA architecture to 7.0 for Tesla V100
set_property(TARGET matrix PROPERTY CUDA_ARCHITECTURES 70)

# Link to Base
target_link_libraries(matrix base)