cmake_minimum_required(VERSION 3.9)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build/cmake/modules/")

include(CheckIPOSupported)

project(octopus)

set(octopus_VERSION_MAJOR 0)
set(octopus_VERSION_MINOR 7)
set(octopus_VERSION_PATCH 4)
set(octopus_VERSION_RELEASE "")

# Generate list of compile commands.  This helps debugging and doesn't have a downside.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Avoid warnings when using GMP_ROOT
cmake_policy(SET CMP0074 NEW)

# Get the current working branch
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH_
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH_
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(GIT_BRANCH "${GIT_BRANCH_}")
set(GIT_COMMIT_HASH "${GIT_COMMIT_HASH_}")

configure_file (
    "${PROJECT_SOURCE_DIR}/src/config/version.h.in"
    "${PROJECT_BINARY_DIR}/generated/version.hpp"
)

include_directories(${CMAKE_BINARY_DIR}/generated)

option(BUILD_SHARED_LIBS "Build the shared library" ON)
set(COMPILER_ARCHITECTURE "native" CACHE STRING "Compiler -march argument")

set(CMAKE_COLOR_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_COMPILER_IS_GNUCXX)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.3")
        message(FATAL_ERROR "GCC version must be at least 6.3!")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.8")
        message(FATAL_ERROR "Clang version must be at least 3.8!")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.8")
        message(FATAL_ERROR "Intel version must be at least 19.0!")
    endif()
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

message("-- Build type: " ${CMAKE_BUILD_TYPE})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${octopus_SOURCE_DIR}/bin" CACHE PATH "Installing to ${octopus_SOURCE_DIR}$" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

message("-- Installation prefix: " ${CMAKE_INSTALL_PREFIX})

if (NOT ${COMPILER_ARCHITECTURE} STREQUAL "native")
    set(CMAKE_SYSTEM_PROCESSOR ${COMPILER_ARCHITECTURE})
endif()
message("-- Target architecture: " ${CMAKE_SYSTEM_PROCESSOR})

# for the main octopus executable
add_subdirectory(lib)
add_subdirectory(src)

# for generating tests
if (BUILD_TESTING)
    add_subdirectory(test)
endif(BUILD_TESTING)
