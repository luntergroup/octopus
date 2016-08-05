set(CONFIG_TEST_SOURCES
    
)

set(CONCEPTS_TEST_SOURCES
    
)

set(BASICS_TEST_SOURCES
    basics/contig_region_tests.cpp
)

set(CONTAINERS_TEST_SOURCES
    
)

set(LOGGING_TEST_SOURCES
    
)

set(IO_TEST_SOURCES
    
)

set(READPIPE_TEST_SOURCES
    
)

set(UTILS_TEST_SOURCES
    
)

set(CORE_TEST_SOURCES
    
)

set(OCTOPUS_TEST_SOURCES
    unit_test_main.cpp
    ${CONFIG_TEST_SOURCES}
    ${CONCEPTS_TEST_SOURCES}
    ${BASICS_TEST_SOURCES}
    ${CONTAINERS_TEST_SOURCES}
    ${LOGGING_TEST_SOURCES}
    ${IO_TEST_SOURCES}
    ${READPIPE_TEST_SOURCES}
    ${UTILS_TEST_SOURCES}
    ${CORE_TEST_SOURCES}
)

add_definitions(-DBOOST_TEST_DYN_LINK)
find_package(Boost 1.58 REQUIRED COMPONENTS unit_test_framework REQUIRED)

include_directories(${Boost_INCLUDE_DIRS} ${octopus_SOURCE_DIR}/lib ${octopus_SOURCE_DIR}/src)

set(TEST_DEPENDENCY_LIBS
    Octopus
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build/cmake/modules/")

include(BoostTestHelpers)
# Add each test
foreach(SRC ${OCTOPUS_TEST_SOURCES})
    add_boost_test(${SRC} ${TEST_DEPENDENCY_LIBS})
endforeach()

# add_executable(test_suite ${OCTOPUS_TEST_SOURCES})
# target_link_libraries(test_suite ${Boost_LIBRARIES})
# install(TARGETS test_suite DESTINATION ${octopus_SOURCE_DIR}/bin)
