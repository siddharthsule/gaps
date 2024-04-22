# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Enable CUDA
enable_language(CUDA) 

# Project name and languages used
project(gaps LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the --expt-extended-lambda flag to the CUDA flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")

# List of subdirectories
set(SUBDIRS base matrix shower observables)

# Include the directories for the headers
foreach(subdir ${SUBDIRS})
  include_directories(${subdir}/include)
  add_subdirectory(${subdir})
endforeach()

# Set the directory for the executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add executable
add_executable(gaps main.cu)

# Link the libraries from the subdirectories
target_link_libraries(gaps base matrix shower observables)

# Set CUDA architecture to 7.0 for Tesla V100
set_property(TARGET gaps PROPERTY CUDA_ARCHITECTURES 70)

# Set CUDA properties
set_target_properties(
    gaps
    PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)