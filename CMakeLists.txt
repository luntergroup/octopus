cmake_minimum_required(VERSION 3.5)

project(octopus)

set(CMAKE_COLOR_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_COMPILER_IS_GNUCXX)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.2")
        message(FATAL_ERROR "Insufficient GCC version (6.2 or greater required)")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.8")
        message(FATAL_ERROR "Insufficient Clang version (3.8 or greater required)")
    endif()
endif()

message("-- Build type is " ${CMAKE_BUILD_TYPE})

# for the main octopus executable
add_subdirectory(lib)
add_subdirectory(src)

# for generating tests
if (BUILD_TESTING)
    add_subdirectory(test)
endif(BUILD_TESTING)
