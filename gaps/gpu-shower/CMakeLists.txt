# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Enable CUDA
enable_language(CUDA)

# Project name and languages used
project(gpu-shower LANGUAGES CXX CUDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the --expt-extended-lambda flag to the CUDA flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda -Xcudafe --diag_suppress=20050")

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/base/include
    ${CMAKE_SOURCE_DIR}/matrix/include
    ${CMAKE_SOURCE_DIR}/shower/include
    ${CMAKE_SOURCE_DIR}/observables/include
)

# Source files
file(GLOB BASE_SOURCES "${CMAKE_SOURCE_DIR}/base/src/*.cu")
file(GLOB MATRIX_SOURCES "${CMAKE_SOURCE_DIR}/matrix/src/*.cu")
file(GLOB SHOWER_SOURCES "${CMAKE_SOURCE_DIR}/shower/src/*.cu")
file(GLOB OBSERVABLES_SOURCES "${CMAKE_SOURCE_DIR}/observables/src/*.cu")
file(GLOB MAIN_SOURCES "${CMAKE_SOURCE_DIR}/main.cu")

# Combine all source files into a single list
set(SOURCES
    ${BASE_SOURCES}
    ${MATRIX_SOURCES}
    ${SHOWER_SOURCES}
    ${OBSERVABLES_SOURCES}
    ${PDF_WRAPPER_SOURCES}
    ${MAIN_SOURCES}
)

# Set all source files to compile with CUDA
set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CUDA)


# Obtain LHAPDF Path
file(READ ${CMAKE_SOURCE_DIR}/../lhapdf-gpu/path.txt LHAPDF_PATH)
string(STRIP ${LHAPDF_PATH} LHAPDF_PATH)

# Add LHAPDF
include_directories(${LHAPDF_PATH}/include)
link_directories(${LHAPDF_PATH}/lib)

# Set the directory for the executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Add executable
add_executable(gpu-shower ${SOURCES})

# Set CUDA architecture to 7.0 for Tesla V100
set_property(TARGET gpu-shower PROPERTY CUDA_ARCHITECTURES 70)

# Set CUDA properties
set_target_properties(
    gpu-shower
    PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Link LHAPDF
target_link_libraries(gpu-shower LHAPDF)