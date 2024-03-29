set(RANGER_SOURCES
    Data.h
    Data.cpp
    DataChar.h
    DataDouble.h
    DataFloat.h
    Forest.h
    Forest.cpp
    ForestClassification.h
    ForestClassification.cpp
    ForestProbability.h
    ForestProbability.cpp
    ForestRegression.h
    ForestRegression.cpp
    ForestSurvival.h
    ForestSurvival.cpp
    globals.h
    Tree.h
    Tree.cpp
    TreeClassification.h
    TreeClassification.cpp
    TreeProbability.h
    TreeProbability.cpp
    TreeRegression.h
    TreeRegression.cpp
    TreeSurvival.h
    TreeSurvival.cpp
    utility.h
    utility.cpp)

set(REQUIRED_BOOST_LIBRARIES
    iostreams
)

find_package (Boost 1.65 REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBRARIES} REQUIRED)

add_library(ranger STATIC ${RANGER_SOURCES})

set(WarningIgnores
    -Wno-unused-parameter
    -Wno-unused-function
    -Wno-missing-braces)

add_compile_options(-Wall -Wextra -Werror ${WarningIgnores})
target_link_libraries (ranger ${Boost_LIBRARIES})
