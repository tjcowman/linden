cmake_minimum_required(VERSION 3.15)

# Must be called before project call
set(CMAKE_GENERATOR_TOOLSET "ClangCL")
set(CMAKE_CXX_COMPILER clang++)

# set the project name
project(linden VERSION 1.0)



# add the executable
add_executable(linden src/main.cpp 
    src/LDForest.cpp 
    src/LDTree2.cpp 
    src/Snp.cpp 
    src/TopSnpList.cpp 
    src/InputParser.cpp "include/Bitwise.h" "include/Types.h" "include/CTable.h"  "include/Graph.h" "src/Graph.cpp")  


set_property(TARGET linden PROPERTY CXX_STANDARD 17)


add_executable(linden-tests test/main.cpp
    "include/Graph.h")




#if(NOT MSCV)
  # target_compile_options( linden PUBLIC -Wall -Wno-c++98-compat-pedantic  -Wno-newline-eof)
#endif()



set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenMP REQUIRED)
target_link_libraries(linden PRIVATE OpenMP::OpenMP_CXX)

find_package(Catch2 3 REQUIRED)
target_link_libraries(linden-tests PRIVATE Catch2::Catch2WithMain)

include_directories(include)

install(TARGETS linden DESTINATION bin)
