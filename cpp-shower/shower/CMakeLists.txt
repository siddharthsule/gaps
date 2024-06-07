cmake_minimum_required(VERSION 3.10)
project(shower)

set(CMAKE_CXX_STANDARD 17)

include_directories(include ../base/include)
file(GLOB SOURCES "src/*.cpp")

add_library(shower ${SOURCES})

# Link to Base
target_link_libraries(shower base)