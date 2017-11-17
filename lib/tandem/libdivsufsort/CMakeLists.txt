set(SOURCES
    config.h
    lfs.h
    divsufsort_private.h
    divsufsort.h
    divsufsort.c
    sssort.c
    trsort.c
    utils.c
)

add_library(libdivsufsort STATIC ${SOURCES})

set_target_properties(libdivsufsort PROPERTIES LINKER_LANGUAGE C)

check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_property(TARGET libdivsufsort PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()