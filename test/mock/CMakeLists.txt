set(MOCK_SOURCES
    mock_reference.hpp
    mock_reference.cpp
)

add_library(Mock ${MOCK_SOURCES})

target_include_directories(Mock PUBLIC ${octopus_SOURCE_DIR}/lib ${octopus_SOURCE_DIR}/src)

target_link_libraries(Mock Octopus)
