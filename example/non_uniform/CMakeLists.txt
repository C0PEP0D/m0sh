cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(non_uniform)

# Build target
add_executable(non_uniform main.cpp)
target_include_directories(non_uniform PRIVATE .)

# Link modules
target_include_directories(non_uniform PRIVATE ../../include)
# Link thirdparty
target_include_directories(non_uniform PRIVATE ../../../../thirdparty/eigen)
