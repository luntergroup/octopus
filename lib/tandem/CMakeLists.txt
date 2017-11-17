add_subdirectory(libdivsufsort)

add_library(tandem STATIC tandem.hpp tandem.cpp)

set(WarningIgnores
    -Wno-unused-parameter
    -Wno-unused-function
    -Wno-missing-braces
    )

add_compile_options(-Wall -Wextra -Werror ${WarningIgnores})
target_link_libraries(tandem libdivsufsort)
set_target_properties(tandem PROPERTIES LINKER_LANGUAGE CXX)

check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_property(TARGET tandem PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
