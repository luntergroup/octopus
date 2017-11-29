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

if(CMAKE_COMPILER_IS_GNUCC)
    target_compile_options(libdivsufsort PRIVATE -fomit-frame-pointer)
endif()
