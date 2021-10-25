cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(unregular)

# Build target
add_executable(unregular main.cpp)
target_include_directories(unregular PRIVATE .)

# Link modules
target_include_directories(unregular PRIVATE ../../include)
# Link thirdparty
target_include_directories(unregular PRIVATE ../../../../thirdparty/eigen)
