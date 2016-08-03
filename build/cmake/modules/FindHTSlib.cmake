# - Try to find htslib 
# Once done, this will define
#
#  htslib_FOUND - system has htslib
#  htslib_INCLUDE_DIRS - the htslib include directories
#  htslib_LIBRARIES - link these to use htslib

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
    $ENV{HTLSIB_ROOT}
    /gsc/pkg/bio/htslib
    /usr
    /usr/local
)

set(_htslib_ver_path "htslib-${htslib_FIND_VERSION}")

# Dependencies
libfind_package(HTSlib ZLIB)

# Include dir
find_path(HTSlib_INCLUDE_DIR
    NAMES ${HTSLIB_ADDITIONAL_HEADERS} sam.h
    PATHS ${HTSLIB_SEARCH_DIRS}
    PATH_SUFFIXES 
        include include/htslib htslib/${_htslib_ver_path}/htslib
    HINTS ENV HTSLIB_ROOT
)

# Finally the library itself
find_library(HTSlib_LIBRARY
    NAMES hts libhts.a hts.a
    PATHS ${HTSlib_INCLUDE_DIR} ${HTSLIB_SEARCH_DIRS}
    NO_DEFAULT_PATH
    PATH_SUFFIXES lib lib64 ${_htslib_ver_path}
    HINTS ENV HTSLIB_ROOT
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(HTSlib_PROCESS_INCLUDES HTSlib_INCLUDE_DIR ZLIB_INCLUDE_DIR)
set(HTSlib_PROCESS_LIBS HTSlib_LIBRARY ZLIB_LIBRARIES)
libfind_process(HTSlib)
message(STATUS "   HTSlib include dirs: ${HTSlib_INCLUDE_DIRS}")
message(STATUS "   HTSlib libraries: ${HTSlib_LIBRARIES}")