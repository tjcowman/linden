cmake_minimum_required(VERSION 3.20)

# Must be called before project call
#set(CMAKE_GENERATOR_TOOLSET "ClangCL")
#set(CMAKE_CXX_COMPILER clang++)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(linden VERSION 1.0)


# add the primary executable
add_executable(linden src/main.cpp 
    src/LDForest.cpp 
    src/LDTree2.cpp 
    src/Snp.cpp 
    src/TopSnpList.cpp 
    src/Bitwise.cpp
    src/SnpSet.cpp
    src/InputParser.cpp "include/Bitwise.h" "include/Types.h" "include/CTable.h"  "include/Graph.h" "src/Graph.cpp" "include/SnpSet.h")  
set_property(TARGET linden PROPERTY CXX_STANDARD 17)


add_executable(linden-serialize src/linden-serialize.cpp
    src/LDForest.cpp 
    src/LDTree2.cpp 
    src/Snp.cpp 
    src/TopSnpList.cpp 
    src/Bitwise.cpp
    src/SnpSet.cpp
    src/InputParser.cpp
    "include/SnpSet.h")
set_property(TARGET linden-serialize PROPERTY CXX_STANDARD 17)

# add the original LDTree version
add_executable(lindenV0 src/main.cpp 
    src/LDForest.cpp 
    src/LDTree.cpp 
    src/Snp.cpp 
    src/TopSnpList.cpp 
    src/Bitwise.cpp
    src/InputParser.cpp "include/Bitwise.h" "include/Types.h" "include/CTable.h"  "include/Graph.h" "src/Graph.cpp" "include/LDTree.h")  
set_property(TARGET lindenV0 PROPERTY CXX_STANDARD 17)
target_compile_definitions(lindenV0 PUBLIC -DoldTree)


add_executable(linden-tests test/main.cpp
    src/LDForest.cpp 
    src/LDTree2.cpp 
    src/Snp.cpp 
    src/TopSnpList.cpp 
    src/Bitwise.cpp
    "include/Graph.h")
set_property(TARGET linden-tests PROPERTY CXX_STANDARD 17)


find_package(OpenMP REQUIRED)


target_link_libraries(linden PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(lindenV0 PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(linden-serialize PRIVATE OpenMP::OpenMP_CXX)


find_package(Catch2 3 REQUIRED)
target_link_libraries(linden-tests PRIVATE Catch2::Catch2WithMain)



include_directories(include)
install(TARGETS linden DESTINATION bin)
