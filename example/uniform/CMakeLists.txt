cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(uniform)

# Build target
add_executable(uniform main.cpp)
target_include_directories(uniform PRIVATE .)

# Link modules
target_include_directories(uniform PRIVATE ../../include)
# Link thirdparty
target_include_directories(uniform PRIVATE ../../../../thirdparty/eigen)
