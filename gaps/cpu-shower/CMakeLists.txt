# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name and languages used
project(cpu-shower LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/base/include
    ${CMAKE_SOURCE_DIR}/matrix/include
    ${CMAKE_SOURCE_DIR}/shower/include
    ${CMAKE_SOURCE_DIR}/observables/include
)

# Optionally enable profiling with -pg flag
option(GPROF "Enable profiling with -pg" OFF)
if(GPROF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

# Source files
file(GLOB BASE_SOURCES "${CMAKE_SOURCE_DIR}/base/src/*.cpp")
file(GLOB MATRIX_SOURCES "${CMAKE_SOURCE_DIR}/matrix/src/*.cpp")
file(GLOB SHOWER_SOURCES "${CMAKE_SOURCE_DIR}/shower/src/*.cpp")
file(GLOB OBSERVABLES_SOURCES "${CMAKE_SOURCE_DIR}/observables/src/*.cpp")
file(GLOB MAIN_SOURCES "${CMAKE_SOURCE_DIR}/main.cpp")

# Combine all source files into a single list
set(SOURCES
    ${BASE_SOURCES}
    ${MATRIX_SOURCES}
    ${SHOWER_SOURCES}
    ${OBSERVABLES_SOURCES}
    ${PDF_WRAPPER_SOURCES}
    ${MAIN_SOURCES}
)

# Obtain LHAPDF Path
file(READ ${CMAKE_SOURCE_DIR}/../lhapdf-gpu/path.txt LHAPDF_PATH)
string(STRIP ${LHAPDF_PATH} LHAPDF_PATH)

# Add LHAPDF
include_directories(${LHAPDF_PATH}/include)
link_directories(${LHAPDF_PATH}/lib)

# Set the directory for the executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add executable
add_executable(cpu-shower ${SOURCES})

# Link LHAPDF
target_link_libraries(cpu-shower LHAPDF)
