cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(sub)

# Build target
add_executable(sub main.cpp)
target_include_directories(sub PRIVATE .)

# Link modules
target_include_directories(sub PRIVATE ../../include)
# Link thirdparty
target_include_directories(sub PRIVATE ../../../../thirdparty/eigen)
