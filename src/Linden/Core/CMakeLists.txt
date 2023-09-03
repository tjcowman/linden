
add_library(Linden-Core STATIC
    Bitwise.h
    Graph.h
    LDForest.cpp
    LDForest.h
    LDTree.cpp
    LDTree.h
    Snp.cpp
    Snp.hpp
    SnpSet.cpp
    SnpSet.h
    TopSnpList.cpp
    TopSnpList.h)
set_property(TARGET Linden-Core PROPERTY CXX_STANDARD 17)
target_link_libraries(Linden-Core PUBLIC
    OpenMP::OpenMP_CXX
    Statistics
    Genetics
    Linden::Load)
target_include_directories(Linden-Core PRIVATE ../LindenLoad)
target_include_directories(Linden-Core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_library(Linden::Core ALIAS Linden-Core)

add_library(Linden-Core-Legacy STATIC
    Bitwise.h
    Graph.h
    LDForest.cpp
    LDForest.h
    Legacy/LDTree.cpp
    Legacy/LDTree.h
    Snp.cpp
    Snp.hpp
    SnpSet.cpp
    SnpSet.h
    TopSnpList.cpp
    TopSnpList.h)
set_property(TARGET Linden-Core-Legacy PROPERTY CXX_STANDARD 17)
target_link_libraries(Linden-Core-Legacy PUBLIC
    OpenMP::OpenMP_CXX
    Statistics
    Genetics
    Linden::Load)
target_include_directories(Linden-Core-Legacy PRIVATE ../LindenLoad)
target_include_directories(Linden-Core-Legacy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_library(Linden::Core::Legacy ALIAS Linden-Core-Legacy)

RegisterTest(Snp.test.cpp)
