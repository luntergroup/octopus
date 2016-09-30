add_subdirectory(libdivsufsort)

add_library(tandem STATIC tandem.hpp tandem.cpp)

set(WarningIgnores
    -Wno-unused-parameter
    -Wno-unused-function
    -Wno-missing-braces
    )

if (CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(tandem PRIVATE -Wall -Wextra -Wextra ${WarningIgnores} -pedantic)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(tandem PRIVATE -Wall -Wextra ${WarningIgnores})
endif()

target_link_libraries(tandem libdivsufsort)

set_target_properties(tandem PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(tandem PROPERTIES INTERPROCEDURAL_OPTIMIZATION True)
