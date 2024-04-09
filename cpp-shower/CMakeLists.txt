# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name and languages used
project(cpp-shower LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# List of subdirectories
set(SUBDIRS base matrix shower observables)

# Include the directories for the headers
foreach(subdir ${SUBDIRS})
  include_directories(${subdir}/include)
  add_subdirectory(${subdir})
endforeach()

# Set the directory for the executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add main.cpp to the executable
add_executable(cpp-shower main.cpp)

# Link the libraries from the subdirectories
target_link_libraries(cpp-shower matrix shower observables)