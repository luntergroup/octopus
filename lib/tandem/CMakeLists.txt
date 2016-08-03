add_subdirectory(libdivsufsort)

add_library(tandem STATIC tandem.hpp tandem.cpp)

target_link_libraries(tandem libdivsufsort)

set_target_properties(tandem PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(tandem PROPERTIES INTERPROCEDURAL_OPTIMIZATION True)
