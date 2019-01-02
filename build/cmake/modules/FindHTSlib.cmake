# - Try to find htslib
# Once done, this will define
#
#  htslib_FOUND - system has htslib
#  htslib_INCLUDE_DIRS - the htslib include directories
#  htslib_LIBRARIES - link these to use htslib
#
# This code was modified from https://github.com/genome/build-common/blob/master/cmake/FindHTSlib.cmake
#

# A simple wrapper to make pkg-config searches a bit easier.
# Works the same as CMake's internal pkg_check_modules but is always quiet.
macro (libfind_pkg_check_modules)
    find_package(PkgConfig QUIET)
    if (PKG_CONFIG_FOUND)
        pkg_check_modules(${ARGN} QUIET)
    endif()
endmacro()

macro (libfind_package PREFIX)
    set (LIBFIND_PACKAGE_ARGS ${ARGN})
    if (${PREFIX}_FIND_QUIETLY)
        set (LIBFIND_PACKAGE_ARGS ${LIBFIND_PACKAGE_ARGS} QUIET)
    endif (${PREFIX}_FIND_QUIETLY)
    if (${PREFIX}_FIND_REQUIRED)
        set (LIBFIND_PACKAGE_ARGS ${LIBFIND_PACKAGE_ARGS} REQUIRED)
    endif (${PREFIX}_FIND_REQUIRED)
    find_package(${LIBFIND_PACKAGE_ARGS})
endmacro (libfind_package)

macro (libfind_process PREFIX)
    # Skip processing if already processed during this run
    if (NOT ${PREFIX}_FOUND)
        # Start with the assumption that the library was found
        set (${PREFIX}_FOUND TRUE)

        # Process all includes and set _FOUND to false if any are missing
        foreach (i ${${PREFIX}_PROCESS_INCLUDES})
            if (${i})
                set (${PREFIX}_INCLUDE_DIRS ${${PREFIX}_INCLUDE_DIRS} ${${i}})
                mark_as_advanced(${i})
            else (${i})
                set (${PREFIX}_FOUND FALSE)
            endif (${i})
        endforeach (i)

        # Process all libraries and set _FOUND to false if any are missing
        foreach (i ${${PREFIX}_PROCESS_LIBS})
            if (${i})
                set (${PREFIX}_LIBRARIES ${${PREFIX}_LIBRARIES} ${${i}})
                mark_as_advanced(${i})
            else (${i})
                set (${PREFIX}_FOUND FALSE)
            endif (${i})
        endforeach (i)

        # Print message and/or exit on fatal error
        if (${PREFIX}_FOUND)
            if (NOT ${PREFIX}_FIND_QUIETLY)
                message (STATUS "Found ${PREFIX} ${${PREFIX}_VERSION}")
            endif (NOT ${PREFIX}_FIND_QUIETLY)
        else (${PREFIX}_FOUND)
            if (${PREFIX}_FIND_REQUIRED)
                foreach (i ${${PREFIX}_PROCESS_INCLUDES} ${${PREFIX}_PROCESS_LIBS})
                    message("${i}=${${i}}")
                endforeach (i)
                message (FATAL_ERROR "Required library ${PREFIX} NOT FOUND.\nInstall the library (dev version) and try again. If the library is already installed, use ccmake to set the missing variables manually.")
            endif (${PREFIX}_FIND_REQUIRED)
        endif (${PREFIX}_FOUND)
    endif (NOT ${PREFIX}_FOUND)
endmacro (libfind_process)

set(HTSLIB_SEARCH_DIRS
    ${HTSLIB_SEARCH_DIRS}
    $ENV{HTSLIB_ROOT}
    ${HTSLIB_ROOT}
    )

if(NOT HTSlib_NO_SYSTEM_PATHS)
    set(HTSLIB_SEARCH_DIRS
        ${HTSLIB_SEARCH_DIRS}
        /usr
        /usr/local
        )
endif()

set(_htslib_ver_path "htslib-${htslib_FIND_VERSION}")

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(HTSLIB_PKGCONF htslib)

# Include dir
find_path(HTSlib_INCLUDE_DIR
    NAMES ${HTSLIB_ADDITIONAL_HEADERS} htslib/sam.h
    PATHS ${HTSLIB_SEARCH_DIRS} ${HTSLIB_PKGCONF_INCLUDE_DIRS}
    PATH_SUFFIXES include htslib/${_htslib_ver_path}
    NO_DEFAULT_PATH
    )

if (HTSlib_USE_STATIC_LIBS)
    # Dependencies
    libfind_package(HTSlib ZLIB)
    libfind_package(HTSlib BZip2)
    libfind_package(HTSlib LibLZMA)
    libfind_package(HTSlib CURL)
    if (NOT APPLE)
        libfind_package(HTSlib OpenSSL)
    endif()
    set(HTSlib_LIBRARY_names libhts.a)
else()
    set(HTSlib_LIBRARY_names libhts.so libhts.so.2 libhts.dylib libhts.2.dylib)
endif()

# Finally the library itself
find_library(HTSlib_LIBRARY
    NAMES ${HTSlib_LIBRARY_names}
    PATHS ${HTSlib_INCLUDE_DIR} ${HTSLIB_SEARCH_DIRS} ${HTSLIB_PKGCONF_LIBRARY_DIRS}
    NO_DEFAULT_PATH
    PATH_SUFFIXES lib lib64 ${_htslib_ver_path}
    )

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(HTSlib_PROCESS_INCLUDES HTSlib_INCLUDE_DIR)
set(HTSlib_PROCESS_LIBS HTSlib_LIBRARY)

if (HTSlib_USE_STATIC_LIBS)
    set(HTSlib_PROCESS_INCLUDES ${HTSlib_PROCESS_INCLUDES}
        ZLIB_INCLUDE_DIR
        BZIP2_INCLUDE_DIR
        LIBLZMA_INCLUDE_DIRS
        CURL_INCLUDE_DIRS)
    set(HTSlib_PROCESS_LIBS ${HTSlib_PROCESS_LIBS}
        ZLIB_LIBRARIES
        BZIP2_LIBRARIES
        LIBLZMA_LIBRARIES
        CURL_LIBRARIES)
    if (NOT APPLE)
        set(HTSlib_PROCESS_INCLUDES ${HTSlib_PROCESS_INCLUDES} OPENSSL_INCLUDE_DIR)
        set(HTSlib_PROCESS_LIBS ${HTSlib_PROCESS_LIBS} OPENSSL_LIBRARIES)
    endif()
endif()

libfind_process(HTSlib)

message(STATUS "   HTSlib include dirs: ${HTSlib_INCLUDE_DIRS}")
message(STATUS "   HTSlib libraries: ${HTSlib_LIBRARIES}")
