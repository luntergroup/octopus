set(DATE_SOURCES
    date.h
    tz.h
    tz.cpp
    tz_private.h
    ios.h)

add_library(date-tz STATIC ${DATE_SOURCES})

set(WarningIgnores
    -Wno-unused-parameter
    -Wno-unused-function
    -Wno-missing-braces)

add_compile_options(-Wall -Wextra -Werror ${WarningIgnores})

# find_package(CURL REQUIRED)
# target_include_directories(date-tz SYSTEM PRIVATE ${CURL_INCLUDE_DIRS})
# target_link_libraries(date-tz PRIVATE ${CURL_LIBRARIES})

