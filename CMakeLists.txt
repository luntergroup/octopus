cmake_minimum_required(VERSION 3.5)

project(octopus)

set(CMAKE_COLOR_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# for the main octopus executable
add_subdirectory(lib)
add_subdirectory(src)

option(DEBUG "Build for debug" OFF)

if(DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else()
    set(CMAKE_BUILD_TYPE Release)
endif()

message("Build type is " ${CMAKE_BUILD_TYPE})

# for generating tests
if (BUILD_TESTING)
    add_subdirectory(test)
endif(BUILD_TESTING)
