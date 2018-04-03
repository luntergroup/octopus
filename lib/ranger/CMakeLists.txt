set(RANGER_SOURCES
    Data.h
    Data.cpp
    DataChar.h
    DataChar.cpp
    DataDouble.h
    DataDouble.cpp
    DataFloat.h
    DataFloat.cpp
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

add_library(ranger STATIC ${RANGER_SOURCES})

set(WarningIgnores
    -Wno-unused-parameter
    -Wno-unused-function
    -Wno-missing-braces)

add_compile_options(-Wall -Wextra -Werror ${WarningIgnores})
