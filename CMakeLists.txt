cmake_minimum_required(VERSION 3.3)

project(octopus)

set(CMAKE_COLOR_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)

# for the main Octopus executable
add_subdirectory(lib)
add_subdirectory(src)

# for generating tests
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif(BUILD_TESTING)
