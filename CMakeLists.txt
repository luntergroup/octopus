cmake_minimum_required (VERSION 3.3)

project (Octopus)

# for the main Octopus executable
add_subdirectory(lib)
add_subdirectory(src)

# for generating tests
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif(BUILD_TESTING)
