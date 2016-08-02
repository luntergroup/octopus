cmake_minimum_required (VERSION 3.3)

project (Octopus)

# Add path for custom CMake modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

set (SOURCE
${SOURCE}
${CMAKE_SOURCE_DIR}/src/main.cpp
)

add_executable(Octopus ${SOURCE})
add_subdirectory(${CMAKE_SOURCE_DIR}/src/hello)
target_link_libraries(Octopus hello)

find_package (Boost 1.61 REQUIRED COMPONENTS system filesystem program_options date_time REQUIRED)
if (Boost_FOUND)
    include_directories (${Boost_INCLUDE_DIR})
    target_link_libraries (Octopus ${Boost_LIBRARIES})
endif (Boost_FOUND)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build/cmake/Modules/")

find_package (HTSlib 1.3.1 REQUIRED)
if (HTSlib_FOUND)
    include_directories (${HTSlib_INCLUDE_DIRS})
    target_link_libraries (Octopus ${HTSlib_LIBRARIES})
endif (HTSlib_FOUND)

target_compile_features (Octopus PRIVATE cxx_range_for)
target_compile_features (Octopus PRIVATE cxx_relaxed_constexpr)
target_compile_features (Octopus PRIVATE cxx_thread_local)
target_compile_features (Octopus PRIVATE cxx_return_type_deduction)

set (CMAKE_BUILD_TYPE Release)
